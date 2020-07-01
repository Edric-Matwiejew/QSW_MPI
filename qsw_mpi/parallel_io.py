import h5py
import numpy as np

def __matricize_index(self, index, offset):

    j = np.int64(np.ceil((1 + index + offset)/self.size[0])) - 1
    i = np.int64(index + offset - j*self.size[0])

    return i, j

def save_result(self, filename, walkname, action = 'a'):

    f = h5py.File(filename + '.h5', action, driver = 'mpio', comm = self.MPI_communicator)

    if self.local_result.ndim == 1:

        dset = f.create_dataset(
                walkname,
                (self.size[0], self.size[0]),
                dtype = np.complex128)

        it, jt = __matricize_index(
                self,
                0,
                self.partition_table[self.rank] - 1)

        kt = 0
        if it != 0:
            dset[it:self.size[0],jt] = self.local_result[0:self.size[0] - it]
            jt += 1
            kt = self.size[0] - it

        ib, jb = __matricize_index(
                self,
                self.rho_v.shape[0] - 1,
                self.partition_table[self.rank] - 1)

        kb = self.local_result.shape[0]
        if ib != self.size[0] - 1:
            kb = kb - ib - 1
            dset[0:ib + 1,jb] = self.local_result[kb:self.local_result.shape[0]]
            jb -= 1

        if jt <= jb:
            dset[:,jt:jb + 1] = np.reshape(
                    self.local_result[kt:kb],
                    (self.size[0], jb - jt + 1), order = "F")

        f.close()

    else:

        dset = f.create_dataset(
                walkname,
                (self.local_result.shape[1], self.size[0], self.size[0]) ,
                dtype = np.complex128)

        it, jt = __matricize_index(
                self,
                0,
                self.partition_table[self.rank] - 1)

        kt = 0
        if it != 0:
            dset[:,it:self.size[0],jt] = self.local_result[0:self.size[0] - it,:].T
            jt += 1
            kt = self.size[0] - it

        ib, jb = __matricize_index(
                self,
                self.rho_v.shape[0] - 1,
                self.partition_table[self.rank] - 1)

        kb = self.local_result.shape[0]
        if ib != self.size[0] - 1:
            kb = kb - ib - 1
            dset[:,0:ib + 1,jb] = self.local_result[kb:self.local_result.shape[0],:].T
            jb -= 1

        if jt <= jb:
            dset[:,:,jt:jb + 1] = np.reshape(
                    self.local_result[kt:kb,:].T,
                    (self.local_result.shape[1],self.size[0],jb-jt+1),
                    order="C")

        f.close()

def save_populations(self, filename, popname, action = 'a'):

    local_populations = self.get_local_populations()

    f = h5py.File(filename + '.h5', action, driver = 'mpio', comm = self.MPI_communicator)

    if self.local_result.ndim == 1:

        save_counts = np.array(self.MPI_communicator.allgather(len(local_populations)), dtype = np.int64)

        bounds = np.zeros(self.flock + 1, dtype = np.int64)
        bounds[0] = 0
        bounds[1:self.flock + 1] = save_counts

        for i in range(self.flock):
            bounds[i + 1] = bounds[i + 1] + bounds[i]

        dset = f.create_dataset(
                popname,
                (self.size[0],),
                dtype = np.complex128)

        dset[bounds[self.rank]:bounds[self.rank + 1]] = local_populations

        f.close()

    else:

        save_counts = np.array(self.MPI_communicator.allgather(local_populations.shape[0]))

        bounds = np.zeros(self.flock + 1, dtype = np.int64)
        bounds[0] = 0
        bounds[1:self.flock + 1] = save_counts

        for i in range(self.flock):
            bounds[i + 1] = bounds[i + 1] + bounds[i]

        dset = f.create_dataset(
                popname,
                (self.size[0], local_populations.shape[1]),
                dtype = np.complex128)

        dset[bounds[self.rank]:bounds[self.rank + 1],::] = local_populations

        f.close()

def save_superoperator(self, filename, operatorname, action = 'a'):

    f = h5py.File(filename + '.h5', action, driver = 'mpio', comm = self.MPI_communicator)

    save_counts = np.array(self.MPI_communicator.allgather(len(self.M_values)))

    rs_save_counts = np.array([self.partition_table[i + 1] - self.partition_table[i] for i in range(self.flock)])
    rs_save_counts[-1] += 1

    bounds = self.__counts_to_offsets(save_counts)
    rs_bounds = self.__counts_to_offsets(rs_save_counts)

    # Save row starts.
    dset = f.create_dataset(
            operatorname + '/indptr',
            (rs_bounds[-1],),
            dtype = np.int64)

    if self.rank + 1 == self.flock:
        dset[rs_bounds[self.rank]:rs_bounds[self.rank + 1]] = self.M_row_starts - 1
    else:
        dset[rs_bounds[self.rank]:rs_bounds[self.rank + 1]] = self.M_row_starts[0:-1] - 1

    # Save column_indices.
    dset = f.create_dataset(
            operatorname + '/indices',
            (bounds[-1],),
            dtype = np.int64)

    dset[bounds[self.rank]:bounds[self.rank + 1]] = self.M_col_indexes - 1

    # Save values.
    dset = f.create_dataset(
            operatorname + '/data',
            (bounds[-1],),
            dtype = np.complex128)

    dset[bounds[self.rank]:bounds[self.rank + 1]] = self.M_values

    f[operatorname].attrs['dimensions'] = (self.M_rows, self.M_rows)

    f.close()


