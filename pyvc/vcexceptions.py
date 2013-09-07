'''
A subclass of the Exception class to handle bad increment input in the daterange generator
'''
class SimFileDoesNotExist(Exception):
    def __init__(self, code):
        self.code = code
    def __str__(self):
        return 'The simulation file {} cannot be found.'.format(self.code)

'''
A subclass of the Exception class to handle bad increment input in the daterange generator
'''
class SimFileNotHDF5(Exception):
    def __init__(self, code):
        self.code = code
    def __str__(self):
        return 'The simulation file {} cannot be opened by h5py. It may not be a valid HDF5 file.'.format(self.code)

class NoSimData(Exception):
    def __init__(self, code):
        self.code = code
    def __str__(self):
        return 'The sim data type {} is invalid. Please use an instance of VCSimData.'.format(type(self.code))

class NoSimFile(Exception):
    def __str__(self):
        return 'The instance of VCSimData has no open sim file. Please use VCSimData.open_file("path/to/file").'

class PlotFormatNotSupported(Exception):
    def __init__(self, code):
        self.code = code
    def __str__(self):
        return 'The plot output format {} is not supported. Please use png, pdf or dat.'.format(self.code)