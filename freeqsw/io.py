import numpy as np
from scipy import sparse as sp
import h5py

class File(object):

    def __init__(self, filename):
        self.File = h5py.File(filename + '.qsw',"a")

    def load_csr(self, group):

        csr = sp.csr_matrix(\
                (np.array(self.File[group]['data'], dtype = np.complex128), \
                np.array(self.File[group]['indices'], dtype =np.int32), \
                np.array(self.File[group]['indptr'], dtype =np.int32)), \
                shape = (self.File[group].attrs['size'], self.File[group].attrs['size']))

        return csr

    def H(self):
        return self.load_csr('H')

    def L(self):
        return self.load_csr('L')

    def sources(self):
        if self.File.attrs['# sources'] > 0:
            return (np.array(self.File['sources']['source sites'], dtype = np.int32), \
                    np.array(self.File['sources']['source rates'], dtype = np.float64))
        else:
            return (None, None)

    def sinks(self):
        if self.File.attrs['# sinks'] > 0:
            return (np.array(self.File['sinks']['sink sites'], dtype = np.int32), \
                    np.array(self.File['sinks']['sink rates'], dtype = np.float64))
        else:
            return (None, None)

    def initial_state(self, state):

        if self.File['initial states gen'][str(name)] is 'FreeQSW':
            raise NameError("FreeQSW generated initial states are not saved to file.")

        return np.array(self.File['initial states'][str(name)], dtype = np.complex128)

    def step(self, name):
        return np.array(self.File['steps'][str(name)], dtype = np.complex128)

    def series(self, name):
        return np.array(self.File['series'][str(name)], dtype = np.complex128)

    def list_steps(self, spacing = 15):

        titles = ["step name", "initial state", "omega", "t"]

        step_names = []
        inital_states = []
        omegas = []
        t = []

        for item in self.File['steps']:
            step_names.append(str(item))
            inital_states.append(str(self.File['steps'][str(item)].attrs['initial state']))
            omegas.append(float(self.File['steps'][str(item)].attrs['omega']))
            t.append(float(self.File['steps'][str(item)].attrs['t']))

        data = [titles] + list(zip(step_names, inital_states, omegas, t))

        for i, d in enumerate(data):
            line = '|'.join(str(x).ljust(spacing) for x in d)
            print(line)
            if i == 0:
                print('-' * len(line))

    def list_series(self, spacing = 15):

        titles = ["series name", "initial state", "omega", "steps", "t1", "t2"]

        series_names = []
        inital_states = []
        omegas = []
        steps = []
        t1 = []
        t2 = []

        for item in self.File['series']:
            series_names.append(str(item))
            inital_states.append(str(self.File['series'][str(item)].attrs['initial state']))
            omegas.append(float(self.File['series'][str(item)].attrs['omega']))
            steps.append(int(self.File['series'][str(item)].attrs['steps']))
            t1.append(float(self.File['series'][str(item)].attrs['t1']))
            t2.append(float(self.File['series'][str(item)].attrs['t2']))

        data = [titles] + list(zip(series_names, inital_states, omegas, steps, t1, t2))

        for i, d in enumerate(data):
            line = '|'.join(str(x).ljust(spacing) for x in d)
            print(line)
            if i == 0:
                print('-' * len(line))

    def system_attributes(self):
        system = {
                "size":int(self.File['H'].attrs['size']),
                "H nnz":int(self.File['H'].attrs['nnz']),
                "L nnz":int(self.File['L'].attrs['nnz']),
                "# sources":int(self.File.attrs["# sources"]),
                "# sinks":int(self.File.attrs["# sinks"])}
        return system

    def close(self):
        self.File.close()


def load_system(filename, MPI_communicator):

    rank = MPI_communicator.Get_rank()

    if rank == 0:
        system = File(str(filename))
        H = system.H()
        L = system.L()
        sources = system.sources()
        sinks = system.sinks()
        size = H.shape[0]
        system.close()
    else:
        size = None
        system = None
        sources = (None, None)
        sinks = (None, None)

    size = MPI_communicator.bcast(size, root = 0)
    sources = MPI_communicator.bcast(sources, root=0)
    sinks = MPI_communicator.bcast(sinks, root=0)

    if rank != 0:
        H = sp.csr_matrix((size,size), dtype = np.complex128)
        L = sp.csr_matrix((size,size), dtype = np.complex128)

    H.data = MPI_communicator.bcast(H.data, root = 0)
    H.indices = MPI_communicator.bcast(H.indices, root = 0)
    H.indptr = MPI_communicator.bcast(H.indptr, root = 0)

    L.data = MPI_communicator.bcast(L.data, root = 0)
    L.indices = MPI_communicator.bcast(L.indices, root = 0)
    L.indptr = MPI_communicator.bcast(L.indptr, root = 0)

    return H, L, sources, sinks


def set_file(
        filename,
        H,
        L,
        source_sites,
        source_rates,
        sink_sites,
        sink_rates,
        action = "a"):

    output = h5py.File(filename + ".qsw", action)

    if 'H' not in output:
        H_out = output.create_group('H')
        H_out.attrs['size'] = H.shape[0]
        H_out.attrs['nnz'] = H.count_nonzero()
        H_out.create_dataset('data', data = H.data)
        H_out.create_dataset('indptr', data = H.indptr)
        H_out.create_dataset('indices', data = H.indices)

    elif (int(H.shape[0])!= int(output['H'].attrs['size'])) \
            or (int(H.count_nonzero()) != int(output['H'].attrs['nnz'])):
        raise NameError('Stored H does not match with the current walk.')

    if 'L' not in output:

        L_out = output.create_group('L')
        L_out.attrs['size'] = L.shape[0]
        L_out.attrs['nnz'] = L.count_nonzero()
        L_out.create_dataset('data', data = L.data)
        L_out.create_dataset('indptr', data = L.indptr)
        L_out.create_dataset('indices', data = L.indices)

    elif (int(L.shape[0]) != int(output['L'].attrs['size'])) \
            or (int(L.count_nonzero()) != int(output['L'].attrs['nnz'])):
        raise NameError('Stored L does not match with the current walk.')

    if 'sources' not in output:

        if source_sites[0] is not -1:
            output.attrs['sources'] = True
            output.attrs['# sources'] = source_sites.shape[0]
            sources = output.create_group('sources')
            sources.create_dataset('source sites', data = source_sites)
            sources.create_dataset('source rates', data = source_rates)
        else:
            output.attrs['# sources'] = 0

    elif source_sites.shape[0] != int(output.attrs['# sources']):
            raise NameError('Stored source sites do not match with the current walk.')

    if 'sinks' not in output:
        if sink_sites[0] is not -1:
            output.attrs['sinks'] = True
            output.attrs['# sinks'] = sink_sites.shape[0]
            sinks = output.create_group('sinks')
            sinks.attrs['sink sites'] = len(sink_sites)
            sinks.create_dataset('sink sites', data = sink_sites)
            sinks.create_dataset('sink rates', data = sink_rates)
        else:
            output.attrs['# sinks'] = 0


    elif source_sites.shape[0] != int(output.attrs['# sinks']):
            raise NameError('Stored sink sites do not match with the current walk.')

    if 'initial states' not in output:
        output.create_group('initial states')
        output['initial states'].attrs['counter'] = 1

    if 'steps' not in output:
        output.create_group('steps')
        output['steps'].attrs['counter'] = 1

    if 'series' not in output:
        output.create_group('series')
        output['series'].attrs['counter'] = 1

    output.close()

