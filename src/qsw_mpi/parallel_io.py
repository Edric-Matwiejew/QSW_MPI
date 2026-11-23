import os
from functools import lru_cache

import h5py
import numpy as np
from mpi4py import MPI


# HDF5 1.x compatible representation of NumPy complex128.
complex_dtype = np.dtype([('r', '<f8'), ('i', '<f8')])


@lru_cache(maxsize=1)
def _native_io():
    try:
        from qsw_mpi.fHDF5 import parallel_io
    except ModuleNotFoundError as exc:
        if exc.name != 'qsw_mpi.fHDF5':
            raise
        return None
    return parallel_io


def _pack_blocks(blocks):
    starts, arrays = zip(*blocks)
    offsets = np.asarray(starts, dtype=np.int64).T
    counts = np.asarray([a.shape for a in arrays], dtype=np.int64).T
    nonempty = [a for a in arrays if a.size]
    if not nonempty:
        buffer = np.empty(1, dtype=np.int8)
    elif len(nonempty) == 1:
        buffer = np.ascontiguousarray(nonempty[0]).reshape(-1).view(np.int8)
    else:
        values = np.empty(sum(a.size for a in nonempty), dtype=nonempty[0].dtype)
        first = 0
        for array in nonempty:
            end = first + array.size
            values[first:end].reshape(array.shape)[:] = array
            first = end
        buffer = values.view(np.int8)
    return offsets, counts, buffer


def _create_datasets(handle, datasets, dimensions):
    for name in datasets:
        if name in handle:
            raise ValueError("Dataset already exists: " + name)
    created = []
    try:
        for name, (shape, dtype, _) in datasets.items():
            handle.create_dataset(name, shape=shape, dtype=dtype)
            created.append(name)
        if dimensions is not None:
            group, shape = dimensions
            handle[group].attrs['dimensions'] = shape
    except Exception:
        for name in created:
            del handle[name]
        raise


def _write_partitions(lib, comm, filename, packed):
    file_id = lib.open_file(os.fsencode(filename), comm.py2f())
    if file_id < 0:
        raise OSError("Could not open HDF5 file for collective output: " + filename)
    try:
        for name, (offsets, counts, buffer) in packed.items():
            status = lib.write_dataset(file_id, name.encode('utf-8'), comm.py2f(),
                                       offsets, counts, buffer)
            if status < 0:
                raise OSError("Could not write HDF5 dataset " + name + " in " + filename)
    finally:
        status = lib.close_file(file_id)
        if comm.allreduce(status < 0, op=MPI.LOR):
            raise OSError("Could not close HDF5 file: " + filename)


def _write(comm, filename, action, datasets, dimensions=None):
    filename = filename + '.h5'
    packed = {}
    error = None
    try:
        lib = _native_io()
        if lib is None:
            raise RuntimeError("Parallel output requires the native HDF5 writer. "
                               "Install parallel HDF5 with Fortran support and run "
                               "'python -m pip install . -Ccmake.define.QSW_PARALLEL_IO=ON' "
                               "from the QSW_MPI checkout.")
        for name, (shape, _, blocks) in datasets.items():
            if 0 not in shape:
                packed[name] = _pack_blocks(blocks)
    except Exception as exc:
        error = exc
    if comm.allreduce(error is not None, op=MPI.LOR):
        raise next(exc for exc in comm.allgather(error) if exc is not None)

    # Close the serial metadata handle before opening the file collectively.
    if comm.rank == 0:
        try:
            with h5py.File(filename, action) as handle:
                _create_datasets(handle, datasets, dimensions)
        except Exception as exc:
            error = exc
    error = comm.bcast(error)
    if error is not None:
        raise error

    try:
        _write_partitions(lib, comm, filename, packed)
    except OSError:
        try:
            if comm.rank == 0:
                with h5py.File(filename, 'r+') as handle:
                    for name in datasets:
                        del handle[name]
        finally:
            comm.Barrier()
        raise


def save_result(self, filename, walkname, action='a'):
    """Saves the density matrix at one or more times to an HDF5 dataset.

    :param self: Walk instance.
    :param filename: File basename.
    :param walkname: Dataset name for the density matrix.
    :param action: File access mode, append ('a') or overwrite ('w').
    """
    n = self.size[0]
    local = self.local_result
    col, row = divmod(int(self.partition_table[self.rank]) - 1, n)
    leading = min(len(local), n - row) if row else 0
    columns, trailing = divmod(len(local) - leading, n)
    first_col = col + int(leading > 0)
    times = local.shape[1:]
    blocks = []
    first = 0
    # A vector partition can start and end partway through a matrix column.
    for r, c, rows, cols in ((row, col, leading, 1), (0, first_col, n, columns),
                             (0, first_col + columns, trailing, 1)):
        end = first + rows * cols
        array = local[first:end].reshape((cols, rows) + times).T
        start = (0,) * len(times) + (r, c)
        blocks.append((start, array))
        first = end

    _write(self.MPI_communicator, filename, action,
           {walkname: (times + (n, n), complex_dtype, blocks)})


def _vertex_populations(self, first, local):
    comm = self.MPI_communicator
    n = len(self.vsets)
    bounds = np.arange(comm.size + 1, dtype=np.int64) * n // comm.size
    sites = np.concatenate(self.vsets).astype(np.intp)
    vertices = np.repeat(np.arange(n), [len(v) for v in self.vsets])
    owned = (sites >= first) & (sites < first + len(local))
    sites = sites[owned] - first
    vertices = vertices[owned]

    # Sum split vertex subspaces on the ranks that will write them.
    steps = local.shape[1] if local.ndim == 2 else 1
    local = local.reshape((len(local), steps))
    result = np.empty((bounds[comm.rank + 1] - bounds[comm.rank], steps),
                      dtype=np.complex128, order='F')
    partial = np.empty(n, dtype=np.complex128)
    counts = np.diff(bounds)
    for step in range(steps):
        partial.fill(0)
        np.add.at(partial, vertices, local[sites, step])
        comm.Reduce_scatter(partial, result[:, step], recvcounts=counts, op=MPI.SUM)
    return bounds[comm.rank], result.reshape((len(result),) + self.local_result.shape[1:])


def save_populations(self, filename, popname, action='a'):
    """Saves the vertex populations at one or more times to an HDF5 dataset.

    :param self: Walk instance.
    :param filename: File basename.
    :param popname: Dataset name for the vertex populations.
    :param action: File access mode, append ('a') or overwrite ('w').
    """
    n = self.size[0]
    first = (int(self.partition_table[self.rank]) - 1 + n) // (n + 1)
    local = self._get_local_populations()
    if getattr(self, 'vsets', None) is not None:
        first, local = _vertex_populations(self, first, local)
        n = len(self.vsets)
    shape = (n,) + local.shape[1:]
    start = (first,) + (0,) * (local.ndim - 1)
    _write(self.MPI_communicator, filename, action,
           {popname: (shape, complex_dtype, [(start, local)])})


def save_superoperator(self, filename, operatorname, action='a'):
    """Saves the superoperator to an HDF5 group as a CSR matrix.

    :param self: Walk instance.
    :param filename: File basename.
    :param operatorname: Group name for the CSR matrix.
    :param action: File access mode, append ('a') or overwrite ('w').
    """
    comm = self.MPI_communicator
    start = int(self.SO_row_starts[0]) - 1
    total = comm.bcast(start + int(self.SO_nnz), root=comm.size - 1)
    first_row = int(self.partition_table[self.rank]) - 1
    count = int(self.SO_local_rows) + int(self.rank == comm.size - 1)
    pointers = self.SO_row_starts[:count].astype(np.int64, copy=False) - 1
    columns = self.SO_col_indexes - 1
    datasets = {
        operatorname + '/indptr': ((self.SO_rows + 1,), np.int64, [((first_row,), pointers)]),
        operatorname + '/indices': ((total,), columns.dtype, [((start,), columns)]),
        operatorname + '/data': ((total,), complex_dtype, [((start,), self.SO_values)])
    }
    _write(comm, filename, action, datasets, (operatorname, (self.SO_rows, self.SO_rows)))
