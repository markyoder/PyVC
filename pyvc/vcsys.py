#!/usr/bin/env python
from pyvc import vcexceptions
from pyvc import vcsimdata

#-------------------------------------------------------------------------------
# The base class for Virtual California simulation objects.
#-------------------------------------------------------------------------------
class VCSys(object):
    def __init__(self, sim_data):
        self.sim_data = sim_data
        
        if type(self.sim_data) != vcsimdata.VCSimData:
            raise vcexceptions.NoSimData(self.sim_data)
        
        if self.sim_data.file is None:
            raise vcexceptions.NoSimFile()
        

