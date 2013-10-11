#!/usr/bin/env python
import tables
import Queue
import multiprocessing
from operator import itemgetter
import time
import itertools
import numpy as np

import math
import calendar
from geographiclib.geodesic import Geodesic

#-------------------------------------------------------------------------------
# A class to perform unit conversions.
#-------------------------------------------------------------------------------
class Converter:
        def __init__(self):
                self.lat0 = 31.5
                self.lon0 = -126.0
                self.earth_radius = self.km_m(6371.0)
    
        #-----------------------------------------------------------------------
        # Given two lat/lon points return the distance between them in meters.
        #-----------------------------------------------------------------------
        def surface_distance(self, lat1, lon1, lat2, lon2):
            conv = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)
            return conv["s12"]
    
        #-----------------------------------------------------------------------
        # Given two lat/lon point return cartesian coords with point 1 at the
        # origin and the given coords at point 2.
        #-----------------------------------------------------------------------
        def latlon_xy_2pt(self, lat1, lon1, lat2, lon2):
            conv = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)
            return (conv["s12"]*math.sin(self.deg_rad(conv["azi1"])), conv["s12"]*math.cos(self.deg_rad(conv["azi1"])))
    
        #-----------------------------------------------------------------------
        # Convert a year given in decimal format (ie. 4.234) and return the
        # year, month, day and optionaly hour, minute, second, fractional second
        #-----------------------------------------------------------------------
        def yearDecimal_YearMonthDay(self,yearDecimal,time=False):
                decimal, year = math.modf(yearDecimal)
                decimal, month = math.modf(12.0*decimal)
                decimal, day = math.modf(calendar.monthrange(1972, int(month + 1))[1]*decimal)
                decimal, hour = math.modf(24.0*decimal)
                decimal, minute = math.modf(60.0*decimal)
                decimal, second = math.modf(60.0*decimal)
                
                if time:
                    return int(year), int(month + 1), int(day + 1), int(hour), int(minute), int(second), decimal
                else:
                    return int(year), int(month + 1), int(day + 1)
    
        #-----------------------------------------------------------------------
        # Convert years to seconds.
        #-----------------------------------------------------------------------
        def year_sec(self, year):
            return 365.0 * 24.0 * 60.0 * 60.0 * year
        
        #-----------------------------------------------------------------------
        # Convert seconds to years.
        #-----------------------------------------------------------------------
        def sec_year(self, sec):
            return (1.0/365.0) * (1.0/24.0) * (1.0/60.0) * (1.0/60.0) * sec
        
        #-----------------------------------------------------------------------
        # Convert kilometers to meters.
        #-----------------------------------------------------------------------
        def km_m(self, km):
            return km * 10.0**(3.0)
        
        #-----------------------------------------------------------------------
        # Convert meters to kilometers.
        #-----------------------------------------------------------------------
        def m_km(self, m):
            return m * 10.0**(-3.0)
        
        #-----------------------------------------------------------------------
        # Convert meters squared to kilometers squared.
        #-----------------------------------------------------------------------
        def msq_kmsq(self, msq):
            return msq * 1000.0**(-2.0)
        
        #-----------------------------------------------------------------------
        # Convert degrees to radians.
        #-----------------------------------------------------------------------
        def deg_rad(self, deg):
            return deg * (math.pi/180.0)
        
        #-----------------------------------------------------------------------
        # Convert radians to degrees.
        #-----------------------------------------------------------------------
        def rad_deg(self, rad):
            return rad * (180.0/math.pi)
        
        #-----------------------------------------------------------------------
        # Set the lat and lon to serve as the zero location.
        #-----------------------------------------------------------------------
        def set_latlon0(self, lat0, lon0):
            self.lat0 = lat0
            self.lon0 = lon0
        
        #-----------------------------------------------------------------------
        # Convert a lat/lon point to cartesian coords relative to lat0/lon0.
        #-----------------------------------------------------------------------
        def latlon_xy(self, lat, lon):
            conv = Geodesic.WGS84.Inverse(self.lat0, self.lon0, lat, lon)
            return (conv["s12"]*math.sin(self.deg_rad(conv["azi1"])), conv["s12"]*math.cos(self.deg_rad(conv["azi1"])))
        
        #-----------------------------------------------------------------------
        # Convert x/y points to lat/lon relative to lat0/lon0.
        #-----------------------------------------------------------------------
        def xy_latlon(self, x, y):
            s12 = (x**2.0 + y**2.0)**(0.5)
            azi1 = self.rad_deg(math.atan2(x,y))
            conv = Geodesic.WGS84.Direct(self.lat0, self.lon0, azi1, s12)
            
            return (conv["lat2"], conv["lon2"])
            
#-------------------------------------------------------------------------------
# parse all of the filters and return a string for the plot title
#-------------------------------------------------------------------------------
def get_plot_label(sim_file, event_range=None, section_filter=None, magnitude_filter=None):
    # always say what file it came from
    label_str = ': from file {}'.format(sim_file)
    # do the event range
    if event_range is not None:
        label_str += ', from {} {}-{}'.format(event_range['type'], event_range['filter'][0], event_range['filter'][1])
    # do the magnitude filter
    if magnitude_filter is not None:
        label_str += ', m {}'.format(magnitude_filter)
    # do the section filter. this is tricky
    if section_filter is not None:
        label_str += ', sections '
        # section_ids are just numbers so they are small, print more of them
        if section_filter['type'] == 'section_id':
            max_sections = 5
        else:
            max_sections = 2
        # if the number of sections is less than our max defined above, print a
        # comma seperated list. if not, just print the first and last.
        if len(section_filter['filter']) < max_sections:
            label_str += ','.join([str(x) for x in section_filter['filter']])
        else:
            label_str += '{}...{}'.format(section_filter['filter'][0], section_filter['filter'][-1])

    return label_str

#-------------------------------------------------------------------------------
# calculate binned averages for scatter plots with the x-axis plotted in logpace
#-------------------------------------------------------------------------------
def calculate_averages(x,y):
    num_bins = math.floor(len(x)/100)
    
    if num_bins < 20:
        num_bins = 20
    elif num_bins > 100:
        num_bins = 100
    
    x = np.array(x)
    y = np.array(y)
    
    if np.min(x) == 0:
        bin_min = 1
    else:
        bin_min = math.floor(math.log(np.min(x),10))
    bin_max = math.ceil(math.log(np.max(x),10))
    
    bins = np.logspace(bin_min,bin_max,num=num_bins)
    inds = np.digitize(x, bins)
    
    binned_data = {}

    for n, i in enumerate(inds):
        try:
            binned_data[i].append(y[n])
        except KeyError:
            binned_data[i] = [y[n]]

    x_ave = []
    y_ave = []
    for k in sorted(binned_data.keys()):
        if k != 0:
            x_ave.append(0.5*(bins[k-1]+bins[k]))
            y_ave.append(sum(binned_data[k])/float(len(binned_data[k])))

    return x_ave, y_ave

#-------------------------------------------------------------------------------
# A class for processing event data in parallel. This is only used if a
# simulation is being analyzed for the first time. The quantities calculated
# here are saved in the simulation file.
#-------------------------------------------------------------------------------
class EventDataProcessor(multiprocessing.Process):
    def __init__(self, sim_file_path, work_queue, result_queue):
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        
        #the path to the sim file
        self.sim_file_path = sim_file_path
        
        super(EventDataProcessor, self).__init__()
    
    def run(self):
    
        #the simulation file
        self.sim_file = tables.open_file(self.sim_file_path, 'r')
        
        #the tables we will need to process
        sweep_table = self.sim_file.get_node('/event_sweep_table')
        event_table = self.sim_file.get_node('/event_table')
        block_info_table = self.sim_file.get_node('/block_info_table')
        
        #these getters speed up the caculation of surface rupture length
        pt1getter = itemgetter(3,4,5)
        pt4getter = itemgetter(18,19,20)
        
        while not self.kill_received:
            # get a task
            try:
                event_range = self.work_queue.get_nowait()
            except Queue.Empty:
                break
            
            # do the processing
            print 'processing events {} - {}'.format(*event_range)
            results = {}
            for evnum in range(*event_range):
                areas = {}
                total_slip = 0.0
                slip_records = 0
                surface_rupture_length = 0.0
                involved_sections = []
                for sweep in sweep_table[event_table[evnum]['start_sweep_rec']:event_table[evnum]['end_sweep_rec']]:
                    eleid = sweep['block_id']
                    if block_info_table[eleid]['m_trace_flag_pt1'] > 0:
                        surface_rupture_length += (sum((x-y)**2.0 for x, y in itertools.izip(pt1getter(block_info_table[eleid]),pt4getter(block_info_table[eleid]))))**0.5
                    
                    involved_sections.append(block_info_table[eleid]['section_id'])
                    areas[eleid] = sweep['area']
                    total_slip += sweep['slip']
                    slip_records += 1
        
                results[evnum] = {'average_slip':total_slip/float(slip_records), 'area':sum(areas.values()), 'surface_rupture_length':surface_rupture_length, 'involved_sections':set(involved_sections)}
            
            self.result_queue.put(results)
        
        self.sim_file.close()

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
        
