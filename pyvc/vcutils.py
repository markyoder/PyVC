#!/usr/bin/env python
import tables

class VCSimData(object):
    def __init__(self, file_path=None):
        self.file = None
        if file_path is not None:
            self.open_file(file_path)
    
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if self.file is not None:
            self.file.close()

    def open_file(self, file_path):
        self.file = tables.open_file(file_path)