#!/usr/bin/env python
from . import VCSys
from . import vcexceptions
from operator import itemgetter

#-------------------------------------------------------------------------------
# A class representing the geometry of a Virtual California simulation.
#-------------------------------------------------------------------------------
class VCGeometry(VCSys):
    def __init__(self,sim_data):
        super(VCGeometry, self).__init__(sim_data)
        
        # the table we need from the simulation data
        self.geometry_data = self.sim_data.file.root.block_info_table
        self.model_extents = self.sim_data.file.root.model_extents
        self.base_lat_lon = self.sim_data.file.root.base_lat_lon
    
    @property
    def base_lat(self):
        return self.base_lat_lon[0]
    
    @property
    def base_lon(self):
        return self.base_lat_lon[1]
    
    @property
    def min_lat(self):
        return self.model_extents.col('min_lat')[0]
    
    @property
    def max_lat(self):
        return self.model_extents.col('max_lat')[0]
    
    @property
    def min_lon(self):
        return self.model_extents.col('min_lon')[0]
    
    @property
    def max_lon(self):
        return self.model_extents.col('max_lon')[0]
    
    @property
    def min_x(self):
        return self.model_extents.col('min_x')[0]
    
    @property
    def max_x(self):
        return self.model_extents.col('max_x')[0]
    
    @property
    def min_y(self):
        return self.model_extents.col('min_y')[0]
    
    @property
    def max_y(self):
        return self.model_extents.col('max_y')[0]
    
    @property
    def min_z(self):
        return self.model_extents.col('min_z')[0]
    
    @property
    def max_z(self):
        return self.model_extents.col('max_z')[0]
    
    def get_fault_traces(self):
        traces = {}
        for block in self.geometry_data.where('(m_trace_flag_pt1 != 0) & (m_trace_flag_pt4 != 0)'):
            try:
                traces[block['section_id']].append( (block['m_x_pt1'], block['m_y_pt1'], block['m_z_pt1']) )
                traces[block['section_id']].append( (block['m_x_pt4'], block['m_y_pt4'], block['m_z_pt4']) )
            except KeyError:
                traces[block['section_id']] = [ (block['m_x_pt1'], block['m_y_pt1'], block['m_z_pt1']) ]
                traces[block['section_id']].append( (block['m_x_pt4'], block['m_y_pt4'], block['m_z_pt4']) )
        return traces
    
    def sections_with_elements(self, elements):
        ele_getter = itemgetter(*elements)
        return set(ele_getter(self.geometry_data.read(field='section_id')))
    
    def events_on_section(self, secid):
        try:
            return self.sim_data.file.root.events_by_section._v_children['section_{}'.format(secid)].read()
        except KeyError:
            raise vcexceptions.BadSectionID(secid)

    def get_section_info(self, section_filter=None, section_id=None):
        if section_filter is not None:
            section_info = {}
            for section in section_filter['filter']:
                # get all of the blocks in the section
                bis = self.geometry_data.read_where('{type} == {value}'.format(type=section_filter['type'], value=section))
                # get the number of blocks along the strike and dip (add 1 cause
                # its zero based)
                bas = max(bis[:]['das_id']) + 1
                bad = max(bis[:]['depth_id']) + 1
                if section_filter['type'] == 'section_id':
                    sn = bis[0]['fault_name']
                    sid = section
                else:
                    sn = section
                    sid = bis[0]['section_id']
                section_info[sid] = {'name':sn, 'blocks_along_strike':bas, 'blocks_along_dip':bad}
            return section_info
        elif section_id is not None:
            return {}
        else:
            section_info = {}
            curr_das = None
            curr_sec = None
            for block in self.geometry_data:
                sec = block['section_id']
                das_id = block['das_id']
                depth_id = block['depth_id']
                if sec != curr_sec:
                    #new sec
                    curr_sec = sec
                    max_bad = 0
                    section_info[curr_sec] = {'name':block['fault_name'], 'blocks_along_strike':0, 'blocks_along_dip':max_bad}
                if das_id != curr_das:
                    #new das
                    section_info[curr_sec]['blocks_along_strike'] += 1
                    curr_das = das_id
                if depth_id > max_bad:
                    max_bad = depth_id
                    section_info[curr_sec]['blocks_along_dip'] = max_bad
            return section_info
    
    def __getitem__(self, bid):
        return self.geometry_data[bid]



            
    