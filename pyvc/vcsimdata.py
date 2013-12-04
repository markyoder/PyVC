#!/usr/bin/env python
import numpy as np

import sys
import math
import Queue
import multiprocessing
import time

import tables

import quakelib

#-------------------------------------------------------------------------------
# A class that manages the connection to a simulation data file. This class
# should be used within a with statement. For example:
#
# with VCSimData() as sim_data:
#     [do stuff with sim_data]
#
# This ensures that the __exit__ method is called and the file is closed when no
# longer needed.
#-------------------------------------------------------------------------------
class VCSimData(object):
    def __init__(self, file_path=None):
        self.file = None
        self.file_path = file_path
        
        self.do_event_area = False
        self.do_event_average_slip = False
        self.do_event_surface_rupture_length = False
        self.do_events_by_section = False
        self.do_das_id = False
        self.do_depth_id = False
        
        if self.file_path is not None:
            self.open_file(self.file_path)
    
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        if self.file is not None:
            self.file.close()
    
    @property
    def filename(self):
        if self.file_path is not None:
            return self.file_path.split('/')[-1]
        else:
            return None

    def open_file(self, file_path):
        if self.file_path is None:
            self.file_path = file_path
        self.file = tables.open_file(self.file_path)

        #-----------------------------------------------------------------------
        # check to see if we need to calculate additional data
        #-----------------------------------------------------------------------
        if 'event_area' not in self.file.root.event_table.colnames:
            self.do_event_area = True

        if 'event_average_slip' not in self.file.root.event_table.colnames:
            self.do_event_average_slip = True
    
        if 'event_surface_rupture_length' not in self.file.root.event_table.colnames:
            self.do_event_surface_rupture_length = True
        
        if 'events_by_section' not in self.file.root._v_groups.keys():
            self.do_events_by_section = True
        
        if self.do_event_area or self.do_event_average_slip or self.do_event_surface_rupture_length or self.do_events_by_section:
            self.calculate_additional_event_data()

        if 'das_id' not in self.file.root.block_info_table.colnames:
            self.do_das_id = True
        
        if 'depth_id' not in self.file.root.block_info_table.colnames:
            self.do_depth_id = True
        
        if self.do_das_id or self.do_depth_id:
            self.calculate_additional_block_data()
        
        if 'model_extents' not in self.file.root._v_children.keys():
            self.calculate_model_extents()
    
    def calculate_model_extents(self):
        #-----------------------------------------------------------------------
        # get info from the original file
        #-----------------------------------------------------------------------
        block_info_table = self.file.root.block_info_table
    
        #-----------------------------------------------------------------------
        # calculate the model extents
        #-----------------------------------------------------------------------
        print 'Calculating model extents'
        
        start_time = time.time()
        
        sys_max_z = -sys.float_info.max
        sys_min_z = sys.float_info.max
        sys_max_x = -sys.float_info.max
        sys_min_x = sys.float_info.max
        sys_max_y = -sys.float_info.max
        sys_min_y = sys.float_info.max
        
        for block in block_info_table:
            min_x = min((block['m_x_pt1'], block['m_x_pt2'], block['m_x_pt3'], block['m_x_pt4']))
            max_x = max((block['m_x_pt1'], block['m_x_pt2'], block['m_x_pt3'], block['m_x_pt4']))
        
            min_y = min((block['m_y_pt1'], block['m_y_pt2'], block['m_y_pt3'], block['m_y_pt4']))
            max_y = max((block['m_y_pt1'], block['m_y_pt2'], block['m_y_pt3'], block['m_y_pt4']))
            
            min_z = min((block['m_z_pt1'], block['m_z_pt2'], block['m_z_pt3'], block['m_z_pt4']))
            max_z = max((block['m_z_pt1'], block['m_z_pt2'], block['m_z_pt3'], block['m_z_pt4']))
            
            if min_x < sys_min_x:
                sys_min_x = min_x
            if max_x > sys_max_x:
                sys_max_x = max_x
                
            if min_y < sys_min_y:
                sys_min_y = min_y
            if max_y > sys_max_y:
                sys_max_y = max_y
            
            if min_z < sys_min_z:
                sys_min_z = min_z
            if max_z > sys_max_z:
                sys_max_z = max_z
        
        base_lat_lon_table = self.file.root.base_lat_lon
        conv = quakelib.Conversion(base_lat_lon_table[0], base_lat_lon_table[1])
        
        ne_corner = conv.convert2LatLon( quakelib.Vec3(sys_max_x, sys_max_y, 0.0) )
        sw_corner = conv.convert2LatLon( quakelib.Vec3(sys_min_x, sys_min_y, 0.0) )
    
        self.file.close()
    
        print 'Done! {} seconds'.format(time.time() - start_time)
        
        print 'Creating new tables'
        table_start_time = time.time()
        
        self.file = tables.open_file(self.file_path, 'a')
        
        desc = {
            'min_x':tables.Float64Col(dflt=0.0),
            'max_x':tables.Float64Col(dflt=0.0),
            'min_y':tables.Float64Col(dflt=0.0),
            'max_y':tables.Float64Col(dflt=0.0),
            'min_z':tables.Float64Col(dflt=0.0),
            'max_z':tables.Float64Col(dflt=0.0),
            'min_lat':tables.Float64Col(dflt=0.0),
            'max_lat':tables.Float64Col(dflt=0.0),
            'min_lon':tables.Float64Col(dflt=0.0),
            'max_lon':tables.Float64Col(dflt=0.0)
        }
        
        model_extents = self.file.create_table('/', 'model_extents', desc, 'Model Extents')
        
        model_extents.row.append()
        model_extents.flush()
        
        model_extents.cols.min_x[0] = sys_min_x
        model_extents.cols.max_x[0] = sys_max_x
        model_extents.cols.min_y[0] = sys_min_y
        model_extents.cols.max_y[0] = sys_max_y
        model_extents.cols.min_z[0] = sys_min_z
        model_extents.cols.max_z[0] = sys_max_z
        model_extents.cols.min_lat[0] = sw_corner.lat()
        model_extents.cols.max_lat[0] = ne_corner.lat()
        model_extents.cols.min_lon[0] = sw_corner.lon()
        model_extents.cols.max_lon[0] = ne_corner.lon()
        
        print 'Done! {} seconds'.format(time.time() - table_start_time)
        
        #-----------------------------------------------------------------------
        # close the file and reopen it with the new table
        #-----------------------------------------------------------------------
        self.file.close()
        self.file = tables.open_file(self.file_path)
        
        print 'Total time {} seconds'.format(time.time() - start_time)
        
    def calculate_additional_block_data(self):
        #-----------------------------------------------------------------------
        # get info from the original file
        #-----------------------------------------------------------------------
        block_info_table = self.file.root.block_info_table
    
        #-----------------------------------------------------------------------
        # get the new block data
        #-----------------------------------------------------------------------
        print 'Calculating new block data'
        
        start_time = time.time()
        
        das_ids = []
        depth_ids = []
        
        curr_das = None
        curr_sec = None
        for block in block_info_table:
            sec = block['section_id']
            das = block['m_das_pt1']
            if sec != curr_sec:
                #new sec
                curr_sec = sec
                das_id = -1
                depth_id = -1
            if das != curr_das:
                #new das
                curr_das = das
                das_id += 1
                depth_id = -1
            depth_id += 1
            das_ids.append(das_id)
            depth_ids.append(depth_id)
    
        print 'Done! {} seconds'.format(time.time() - start_time)
        
        print 'Creating new tables'
        table_start_time = time.time()
        #-----------------------------------------------------------------------
        # create the new block_info_table
        #-----------------------------------------------------------------------
        #close the current file
        self.file.close()
        #reopen the file with the append flag set
        self.file = tables.open_file(self.file_path, 'a')
        block_info_table = self.file.root.block_info_table
        
        # get a description of table in dictionary format
        desc_orig = block_info_table.description._v_colObjects
        desc_new = desc_orig.copy()
        
        # add columns to description
        if self.do_das_id:
            desc_new['das_id'] = tables.UIntCol(dflt=0)
        
        if self.do_depth_id:
            desc_new['depth_id'] = tables.UIntCol(dflt=0)
        
        # create a new table with the new description
        block_info_table_new = self.file.create_table('/', 'tmp', desc_new, 'Block Info Table')
        
        # copy the user attributes
        block_info_table.attrs._f_copy(block_info_table_new)
        
        # fill the rows of new table with default values
        for i in xrange(block_info_table.nrows):
            block_info_table_new.row.append()
        
        # flush the rows to disk
        block_info_table_new.flush()
        
        # copy the columns of source table to destination
        for col in desc_orig:
            getattr(block_info_table_new.cols, col)[:] = getattr(block_info_table.cols, col)[:]

        # fill the new columns
        if self.do_das_id:
            block_info_table_new.cols.das_id[:] = das_ids
        
        if self.do_depth_id:
            block_info_table_new.cols.depth_id[:] = depth_ids
        	
        # remove the original table
        block_info_table.remove()
        
        # move table2 to table
        block_info_table_new.move('/','block_info_table')
        
        print 'Done! {} seconds'.format(time.time() - table_start_time)
        
        #-----------------------------------------------------------------------
        # close the file and reopen it with the new tables
        #-----------------------------------------------------------------------
        self.file.close()
        self.file = tables.open_file(self.file_path)
    
        print 'Total time {} seconds'.format(time.time() - start_time)
    
    def calculate_additional_event_data(self):
        #-----------------------------------------------------------------------
        # get info from the original file
        #-----------------------------------------------------------------------
        total_events = self.file.root.event_table.nrows
        #close the current file
        self.file.close()
        
        #-----------------------------------------------------------------------
        # get the new event data
        #-----------------------------------------------------------------------
        print 'Calculating new event data'
        
        start_time = time.time()
        num_processes = multiprocessing.cpu_count()
        
        # break the work up
        seg = int(round(float(total_events)/float(num_processes)))
        work_queue = multiprocessing.Queue()
        for i in range(num_processes):
            if i == num_processes - 1:
                end_index = total_events
            else:
                end_index = seg*int(i + 1)
            work_queue.put((int(i) * seg , end_index))

        # create a queue to pass to workers to store the results
        result_queue = multiprocessing.Queue()

        # spawn workers
        for i in range(num_processes):
            worker = EventDataProcessor(self.file_path, work_queue, result_queue)
            worker.start()

        # collect the results off the queue
        results = {}
        for i in range(num_processes):
            results = dict(results, **result_queue.get())
        data_process_results_sorted = [results[key] for key in sorted(results.keys())]
        
        print 'Done! {} seconds'.format(time.time() - start_time)
        
        print 'Creating new tables'
        table_start_time = time.time()
        #-----------------------------------------------------------------------
        # create the new event_table
        #-----------------------------------------------------------------------
        self.file = tables.open_file(self.file_path, 'a')
        event_table = self.file.root.event_table
        
        # get a description of table in dictionary format
        desc_orig = event_table.description._v_colObjects
        desc_new = desc_orig.copy()
        
        # add columns to description
        if self.do_event_area:
            desc_new['event_area'] = tables.Float64Col(dflt=0.0)
        
        if self.do_event_average_slip:
            desc_new['event_average_slip'] = tables.Float64Col(dflt=0.0)
        
        if self.do_event_surface_rupture_length:
            desc_new['event_surface_rupture_length'] = tables.Float64Col(dflt=0.0)
        
        # create a new table with the new description
        event_table_new = self.file.create_table('/', 'tmp', desc_new, 'Event Table')
        
        # copy the user attributes
        event_table.attrs._f_copy(event_table_new)
        
        # fill the rows of new table with default values
        for i in xrange(event_table.nrows):
            event_table_new.row.append()
        
        # flush the rows to disk
        event_table_new.flush()
        
        # copy the columns of source table to destination
        for col in desc_orig:
            getattr(event_table_new.cols, col)[:] = getattr(event_table.cols, col)[:]

        # fill the new columns
        if self.do_event_area:
            event_table_new.cols.event_area[:] = [ x['area'] for x in data_process_results_sorted ]
        
        if self.do_event_average_slip:
            event_table_new.cols.event_average_slip[:] = [ x['average_slip'] for x in data_process_results_sorted ]

        if self.do_event_surface_rupture_length:
            event_table_new.cols.event_surface_rupture_length[:] = [ x['surface_rupture_length'] for x in data_process_results_sorted ]
        	
        # remove the original table
        event_table.remove()
        
        # move table2 to table
        event_table_new.move('/','event_table')
        
        #-----------------------------------------------------------------------
        # create a new group to store the event ids for each section
        #-----------------------------------------------------------------------
        if self.do_events_by_section:
            events_by_section_group = self.file.create_group("/", 'events_by_section', 'Events on each section')
            
            # dict to store the events on each section
            events_by_section = {}
            
            # we have the sections involved in each event. we need to transpose
            # this to events on each section
            for evid, event_data in enumerate(data_process_results_sorted):
                for secid in event_data['involved_sections']:
                    try:
                        events_by_section[secid].append(evid)
                    except KeyError:
                        events_by_section[secid] = [evid]
            
            # create an array for each result
            for secid, evids in events_by_section.iteritems():
                self.file.create_array(events_by_section_group, 'section_{}'.format(secid), np.array(evids), 'Events on section {}'.format(secid))
            
        print 'Done! {} seconds'.format(time.time() - table_start_time)
        
        #-----------------------------------------------------------------------
        # close the file and reopen it with the new tables
        #-----------------------------------------------------------------------
        self.file.close()
        self.file = tables.open_file(self.file_path)
    
        print 'Total time {} seconds'.format(time.time() - start_time)