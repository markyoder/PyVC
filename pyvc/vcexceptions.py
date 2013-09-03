'''
A subclass of the Exception class to handle bad increment input in the daterange generator
'''
class SimFileDoesNotExist(Exception):
    def __init__(self, code):
        self.code = code
    def __str__(self):
        return 'The simulation file %s cannot be found.'%self.code

'''
A subclass of the Exception class to handle bad increment input in the daterange generator
'''
class SimFileNotHDF5(Exception):
    def __init__(self, code):
        self.code = code
    def __str__(self):
        return 'The simulation file %s cannot be opened by h5py. It may not be a valid HDF5 file.'%self.code