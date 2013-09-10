#!/usr/bin/env python
from . import VCSys

#-------------------------------------------------------------------------------
# A class representing the geometry of a Virtual California simulation.
#-------------------------------------------------------------------------------
class VCGeometry(VCSys):
    def __init__(self,sim_data):
        super(VCGeometry, self).__init__(sim_data)
        
        # the table we need from the simulation data
        self.geometry_data = self.sim_data.file.root.block_info_table
            
    