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

class Converter:
        def __init__(self):
                self.lat0 = 31.5
                self.lon0 = -126.0
                self.earth_radius = self.km_m(6371.0)
        
        def distanceOnEarthsSurface(self, lat1, lon1, lat2, lon2):
            conv = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)
            return conv["s12"]
        
        def latlon_xy_2pt(self, lat1, lon1, lat2, lon2):
            conv = Geodesic.WGS84.Inverse(lat1, lon1, lat2, lon2)
            return (conv["s12"]*math.sin(self.deg_rad(conv["azi1"])), conv["s12"]*math.cos(self.deg_rad(conv["azi1"])))
        
        def yearDecimalToYearMonthDay(self,yearDecimal,time=False):
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
        
        def year_sec(self, year):
            return 365.0 * 24.0 * 60.0 * 60.0 * year
        
        def sec_year(self, sec):
            return (1.0/365.0) * (1.0/24.0) * (1.0/60.0) * (1.0/60.0) * sec
            
        def km_m(self, km):
            return km * 10.0**(3.0)
        
        def m_km(self, m):
            return m * 10.0**(-3.0)

        def msq_kmsq(self, msq):
            return msq * 1000.0**(-2.0)
        
        def deg_rad(self, deg):
            return deg * (math.pi/180.0)
        
        def rad_deg(self, rad):
            return rad * (180.0/math.pi)
    
        def setLatLon0(self, lat0, lon0):
            self.lat0 = lat0
            self.lon0 = lon0

        def latlon_xy(self, lat, lon):
            conv = Geodesic.WGS84.Inverse(self.lat0, self.lon0, lat, lon)
            return (conv["s12"]*math.sin(self.deg_rad(conv["azi1"])), conv["s12"]*math.cos(self.deg_rad(conv["azi1"])))
            
        def latlon_xy_old(self, lat, lon):
            x = self.arclen((self.lat0, self.lon0), (self.lat0, lon))
            y = self.arclen((self.lat0,0),(lat,0))
            
            if (lon < self.lon0):
                x *= -1
            
            if (lat < self.lat0):
                y *= -1
            
            return (x,y)

        def xy_latlon(self, x, y):
            s12 = (x**2.0 + y**2.0)**(0.5)
            azi1 = self.rad_deg(math.atan2(x,y))
            conv = Geodesic.WGS84.Direct(self.lat0, self.lon0, azi1, s12)
            
            return (conv["lat2"], conv["lon2"])
            
        def xy_latlon_old(self, x, y):
            
            new_lat = self.rad_deg(y/self.earth_radius)+self.lat0;

            new_lon = 2.0 * self.rad_deg(math.asin(math.sin(x/(2.0 * self.earth_radius))/math.cos(self.deg_rad(new_lat))))+self.lon0;

            return new_lat,new_lon

        def arclen(self, pt1, pt2):
            # pts are always (lat,lon)

            dlon = self.deg_rad(pt2[1]-pt1[1])
            dlat = self.deg_rad(pt2[0]-pt1[0])
            lat1 = self.deg_rad(pt1[0])
            lat2 = self.deg_rad(pt2[0])

            a = math.sin(dlat/2.0)**2.0 + math.cos(lat1)*math.cos(lat2)*( math.sin(dlon/2.0)**2.0 )
            c = 2.0*math.atan2(math.sqrt(a), math.sqrt(1.0-a))
            d = self.earth_radius*c

            return d

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
                    #trace_sum = block_info_table[eleid]['m_trace_flag_pt1'] + block_info_table[eleid]['m_trace_flag_pt2'] + block_info_table[eleid]['m_trace_flag_pt3'] + block_info_table[eleid]['m_trace_flag_pt4']
                    
                    if block_info_table[eleid]['m_trace_flag_pt1'] > 0:
                        #pt1 = (block_info_table[eleid]['m_x_pt1'], block_info_table[eleid]['m_y_pt1'], block_info_table[eleid]['m_z_pt1'])
                        #pt4 = (block_info_table[eleid]['m_x_pt4'], block_info_table[eleid]['m_y_pt4'], block_info_table[eleid]['m_z_pt4'])
                        surface_rupture_length += (sum((x-y)**2.0 for x, y in itertools.izip(pt1getter(block_info_table[eleid]),pt4getter(block_info_table[eleid]))))**0.5
                    
                    involved_sections.append(block_info_table[eleid]['section_id'])
                    areas[eleid] = sweep['area']
                    total_slip += sweep['slip']
                    slip_records += 1
                #print areas.keys()
                results[evnum] = {'average_slip':total_slip/float(slip_records), 'area':sum(areas.values()), 'surface_rupture_length':surface_rupture_length, 'involved_sections':set(involved_sections)}
            
                #event_elements[i] = set(events.get_event_elements(i))
            
           
            self.result_queue.put(results)
        
        self.sim_file.close()


class VCSimData(object):
    def __init__(self, file_path=None):
        self.file = None
        self.file_path = file_path
        
        self.do_event_area = False
        self.do_event_average_slip = False
        self.do_event_surface_rupture_length = False
        self.do_events_by_section = False
        
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
        
        if 'event_area' not in self.file.root.event_table.colnames:
            self.do_event_area = True

        if 'event_average_slip' not in self.file.root.event_table.colnames:
            self.do_event_average_slip = True
    
        if 'event_surface_rupture_length' not in self.file.root.event_table.colnames:
            self.do_event_surface_rupture_length = True
        
        if 'events_by_section' not in self.file.root._v_groups.keys():
            self.do_events_by_section = True
        
        if self.do_event_area or self.do_event_average_slip or self.do_event_surface_rupture_length or self.do_events_by_section:
            self.calculate_additional_data()

    def calculate_additional_data(self):
        #-----------------------------------------------------------------------
        # get info from the original file
        #-----------------------------------------------------------------------
        total_events = self.file.root.event_table.nrows
        #close the current file
        self.file.close()
        
        #-----------------------------------------------------------------------
        # get the new data
        #-----------------------------------------------------------------------
        print 'Calculating new data'
        
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
        
