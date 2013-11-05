#!/usr/bin/env python
from pyvc import *
from pyvc import vcutils
from pyvc import vcexceptions
import matplotlib.pyplot as mplt
import matplotlib.font_manager as mfont
import matplotlib.colors as mcolor
import matplotlib.colorbar as mcolorbar
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import numpy as np
import math
import multiprocessing
import Queue
import cPickle
import networkx as nx
from operator import itemgetter
from mpl_toolkits.basemap import Basemap, maskoceans, interp
import quakelib
import time
from PIL import Image

class DisplacementGridProcessor(multiprocessing.Process):
    def __init__(self, work_queue, result_queue, field_1d, event_element_data, event_element_slips, lat_size, lon_size, cutoff, event_center, event_radius):#, min_lat, min_lon, max_lat, max_lon):
        
        # base class initialization
        multiprocessing.Process.__init__(self)
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        
        self.field_1d = field_1d
        self.event_element_data = event_element_data
        self.event_element_slips = event_element_slips
        self.lat_size = lat_size
        self.lon_size = lon_size
        self.cutoff = cutoff
        self.event_center = event_center
        self.event_radius = event_radius
    
        #self.counter = counter
        #self.total_tasks = total_tasks
    
    def run(self):
        while not self.kill_received:
            # get a task
            try:
                start, end = self.work_queue.get_nowait()
            except Queue.Empty:
                break
            
            print 'processing elements {} - {}'.format(start, end)
            
            # empty arrays to store the results
            dX = np.empty((self.lat_size, self.lon_size))
            dY = np.empty((self.lat_size, self.lon_size))
            dZ = np.empty((self.lat_size, self.lon_size))
            
            # create a element list
            elements = quakelib.EventElementList()
            
            # create elements and add them to the element list
            for element in self.event_element_data[start:end]:
                #print element
                ele = quakelib.EventElement4()
                ele.set_rake(element['rake_rad'])
                ele.set_slip(self.event_element_slips[element['block_id']])
                ele.set_vert(0, element['m_x_pt1'], element['m_y_pt1'], element['m_z_pt1'])
                ele.set_vert(1, element['m_x_pt2'], element['m_y_pt2'], element['m_z_pt2'])
                ele.set_vert(2, element['m_x_pt3'], element['m_y_pt3'], element['m_z_pt3'])
                ele.set_vert(3, element['m_x_pt4'], element['m_y_pt4'], element['m_z_pt4'])
                elements.append(ele)
            
            # create an event
            event = quakelib.Event()
            
            # add the elements to the event
            event.add_elements(elements)
            
            # set the event center and radius
            event.set_event_center(self.event_center)
            event.set_event_radius(self.event_radius)
            
            #calculate the displacements
            lame_lambda = 3.2e10
            lame_mu = 3.0e10
            
            if self.cutoff is None:
                disp_1d = event.event_displacements(self.field_1d, lame_lambda, lame_lambda)
            else:
                disp_1d = event.event_displacements(self.field_1d, lame_lambda, lame_lambda, self.cutoff)
            disp = np.array(disp_1d).reshape((self.lat_size,self.lon_size))
        
            it = np.nditer(dX, flags=['multi_index'])
            while not it.finished:
                dX[it.multi_index] = disp[it.multi_index][0]
                dY[it.multi_index] = disp[it.multi_index][1]
                dZ[it.multi_index] = disp[it.multi_index][2]
                it.iternext()
            
            # store the result
            processed_displacements = {}
            processed_displacements['dX'] = dX
            processed_displacements['dY'] = dY
            processed_displacements['dZ'] = dZ
            self.result_queue.put(processed_displacements)

#-------------------------------------------------------------------------------
# A class to handle the plotting of event displacements
#-------------------------------------------------------------------------------
class VCDisplacementMapPlotter(object):
    #---------------------------------------------------------------------------
    # If outut_file is none returns an instance for further plotting
    #---------------------------------------------------------------------------
    def __init__(self, min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, fault_traces, padding=0.01, map_res='i', map_proj='cyl'):
        self.look_azimuth = None
        self.look_elevation = None
        self.wavelength = 0.03
        
        # Define how the cutoff value scales if it is not explitly set
        self.cutoff_min_size = 20.0
        self.cutoff_min = 46.5
        self.cutoff_p2_size = 65.0
        self.cutoff_p2 = 90.0
        #self.output_file = None
        #self.min_lon = min_lon
        #self.max_lon = max_lon
        #self.min_lat = min_lat
        #self.max_lat = max_lat
        #self.padding = pading
        
        # These are constrained this way so we can plot on 1024x780 for the
        # animations
        max_map_width = 690.0
        max_map_height = 658.0
        
        self.convert = quakelib.Conversion(base_lat, base_lon)
        
        #-----------------------------------------------------------------------
        # DisplacementmMap configuration
        #-----------------------------------------------------------------------
        # values for the fringes map are denoted by a {value}_f
        self.dmc = {
            'font':               mfont.FontProperties(family='Arial', style='normal', variant='normal', weight='normal'),
            'font_bold':          mfont.FontProperties(family='Arial', style='normal', variant='normal', weight='bold'),
            'cmap':               mplt.get_cmap('YlOrRd'),
            'cmap_f':             mplt.get_cmap('jet'),
        #water
            'water_color':          '#4eacf4',
            'water_color_f':        '#4eacf4',
        #map boundaries
            'boundary_color':       '#000000',
            'boundary_color_f':     '#ffffff',
            'boundary_width':       1.0,
            'coastline_color':      '#000000',
            'coastline_color_f':    '#ffffff',
            'coastline_width':      1.0,
            'country_color':        '#000000',
            'country_color_f':      '#ffffff',
            'country_width':        1.0,
            'state_color':          '#000000',
            'state_color_f':        '#ffffff',
            'state_width':          1.0,
        #rivers
            'river_width':          0.25,
        #faults
            'fault_color':          '#000000',
            'fault_color_f':        '#ffffff',
            'event_fault_color':    '#ff0000',
            'event_fault_color_f':  '#ffffff',
            'fault_width':          0.5,
        #lat lon grid
            'grid_color':           '#000000',
            'grid_color_f':         '#ffffff',
            'grid_width':           0.0,
            'num_grid_lines':       5,
        #map props
            'map_resolution':       map_res,
            'map_projection':       map_proj,
            'plot_resolution':      72.0,
            'map_tick_color':       '#000000',
            'map_tick_color_f':     '#000000',
            'map_frame_color':      '#000000',
            'map_frame_color_f':    '#000000',
            'map_frame_width':      1,
            'map_fontsize':         12,
            'arrow_inset':          10.0,
            'arrow_fontsize':       9.0,
            'cb_fontsize':          10.0,
            'cb_fontcolor':         '#000000',
            'cb_fontcolor_f':       '#000000',
            'cb_height':            20.0,
            'cb_margin_t':          10.0
        }
        
        lon_range = max_lon - min_lon
        lat_range = max_lat - min_lat
        max_range = max((lon_range, lat_range))
        self.padded_min_lon = min_lon - lon_range*padding
        self.padded_min_lat = min_lat - lat_range*padding
        self.padded_max_lon = max_lon + lon_range*padding
        self.padded_max_lat = max_lat + lat_range*padding
        
        #-----------------------------------------------------------------------
        # m1, fig1 is the oceans and the continents. This will lie behind the
        # masked data image.
        #-----------------------------------------------------------------------
        self.m1 = Basemap(
            llcrnrlon=self.padded_min_lon,
            llcrnrlat=self.padded_min_lat,
            urcrnrlon=self.padded_max_lon,
            urcrnrlat=self.padded_max_lat,
            lat_0=(self.padded_max_lat+self.padded_min_lat)/2.0,
            lon_0=(self.padded_max_lon+self.padded_min_lon)/2.0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        #-----------------------------------------------------------------------
        # m2, fig2 is the plotted deformation data.
        #-----------------------------------------------------------------------
        self.m2 = Basemap(
            llcrnrlon=self.padded_min_lon,
            llcrnrlat=self.padded_min_lat,
            urcrnrlon=self.padded_max_lon,
            urcrnrlat=self.padded_max_lat,
            lat_0=(self.padded_max_lat+self.padded_min_lat)/2.0,
            lon_0=(self.padded_max_lon+self.padded_min_lon)/2.0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        #-----------------------------------------------------------------------
        # m3, fig3 is the ocean land mask.
        #-----------------------------------------------------------------------
        self.m3 = Basemap(
            llcrnrlon=self.padded_min_lon,
            llcrnrlat=self.padded_min_lat,
            urcrnrlon=self.padded_max_lon,
            urcrnrlat=self.padded_max_lat,
            lat_0=(self.padded_max_lat+self.padded_min_lat)/2.0,
            lon_0=(self.padded_max_lon+self.padded_min_lon)/2.0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        
        #-----------------------------------------------------------------------
        # Calculate the map grid
        #-----------------------------------------------------------------------
        # aspect is height/width
        
        if self.m1.aspect > 1.0:
            map_height = max_map_height
            map_width = max_map_height/self.m1.aspect
        else:
            map_width = max_map_width
            map_height = max_map_width*self.m1.aspect
        
        #'''
        self.lons_1d = np.linspace(self.padded_min_lon,self.padded_max_lon,int(map_width))
        self.lats_1d = np.linspace(self.padded_min_lat,self.padded_max_lat,int(map_height))
        
        _lons_1d = quakelib.FloatList()
        _lats_1d = quakelib.FloatList()
        
        for lon in self.lons_1d:
            _lons_1d.append(lon)
        
        for lat in self.lats_1d:
            _lats_1d.append(lat)
        
        self.field_1d = self.convert.convertArray2xyz(_lats_1d,_lons_1d)
        #'''
        #self.lats_1d,self.lons_1d,self.field_1d = cPickle.load(open('local/test_grid.pkl','rb'))
        
        self.fault_traces = {}
        for secid in fault_traces.iterkeys():
             self.fault_traces[secid] = zip(*[(lambda y: (y.lat(),y.lon()))(self.convert.convert2LatLon(quakelib.Vec3(x[0], x[1], x[2]))) for x in fault_traces[secid]])
    
    #---------------------------------------------------------------------------
    # Sets up the displacement calculation and then passes it to the
    # DisplacementGridProcessor class.
    #---------------------------------------------------------------------------
    def calculate_displacements(self, event_element_data, event_element_slips, cutoff=None):
    
        #-----------------------------------------------------------------------
        # We need to calculate the event center and radius
        #-----------------------------------------------------------------------
        # create a element list
        elements = quakelib.EventElementList()
        
        # create elements and add them to the element list
        for element in event_element_data:
            ele = quakelib.EventElement4()
            ele.set_vert(0, element['m_x_pt1'], element['m_y_pt1'], element['m_z_pt1'])
            ele.set_vert(1, element['m_x_pt2'], element['m_y_pt2'], element['m_z_pt2'])
            ele.set_vert(2, element['m_x_pt3'], element['m_y_pt3'], element['m_z_pt3'])
            ele.set_vert(3, element['m_x_pt4'], element['m_y_pt4'], element['m_z_pt4'])
            elements.append(ele)
        
        # create an event
        event = quakelib.Event()
        
        # add the elements to the event
        event.add_elements(elements)
        
        # find the event radius and center
        event_center = event.event_center()
        event_radius = event.event_radius()
        
        #print event_center, event_radius
        
        #self.event_center = event_center
        #self.event_radius = event_radius
        
        #-----------------------------------------------------------------------
        # If the cutoff is none (ie not explicitly set) calculate the cutoff for
        # this event.
        #-----------------------------------------------------------------------
        event_size = float(len(event_element_slips))
        if cutoff is None:
            if  event_size >= self.cutoff_min_size:
                cutoff = vcutils.linear_interp(
                    event_size,
                    self.cutoff_min_size,
                    self.cutoff_p2_size,
                    self.cutoff_min,
                    self.cutoff_p2
                    )
            else:
                cutoff = self.cutoff_min
        
        print cutoff
        #-----------------------------------------------------------------------
        # Now break up the work and start the workers
        #-----------------------------------------------------------------------

        # How many seperate CPUs do we have
        num_processes = multiprocessing.cpu_count()
        
        # Figure out how many segments we will break the task up into
        seg = int(round(event_size/float(num_processes)))
        if seg < 1:
            seg = 1
        
        # Break up the job.
        segmented_elements_indexes = []
        for i in range(num_processes):
            if i == num_processes - 1:
                end_index = len(event_element_slips)
            else:
                end_index = seg*int(i + 1)
            start_index = int(i) * seg
            if start_index != end_index:
                segmented_elements_indexes.append((start_index, end_index))
        
        # Add all of the jobs to a work queue
        work_queue = multiprocessing.Queue()
        for job in segmented_elements_indexes:
            work_queue.put(job)
        
        # Create a queue to pass to workers to store the results
        result_queue = multiprocessing.Queue()
        
        # Spawn workers
        for i in range(len(segmented_elements_indexes)):
            worker = DisplacementGridProcessor(work_queue,
                result_queue,
                self.field_1d,
                event_element_data,
                event_element_slips,
                self.lats_1d.size,
                self.lons_1d.size,
                cutoff,
                event_center,
                event_radius
            )
            worker.start()
        
        # Collect the results off the queue
        results = []
        for i in range(len(segmented_elements_indexes)):
            results.append(result_queue.get())
        
        # Combine the results
        self.dX = None
        self.dY = None
        self.dZ = None
        for result_num, result in enumerate(results):
            if self.dX is None:
                self.dX = result['dX']
            else:
                self.dX += result['dX']
            
            if self.dY is None:
                self.dY = result['dY']
            else:
                self.dY += result['dY']
                
            if self.dZ is None:
                self.dZ = result['dZ']
            else:
                self.dZ += result['dZ']

    #---------------------------------------------------------------------------
    # Calculates the look angles based on a set of elements. This will set the
    # angles to be looking along the average strike of the fault.
    #---------------------------------------------------------------------------
    def calculate_look_angles(self, element_data):
        strikes = []
        rakes = []
        
        for element in element_data:
            ele = quakelib.Element4()
            ele.set_vert(0, element['m_x_pt1'], element['m_y_pt1'], element['m_z_pt1'])
            ele.set_vert(1, element['m_x_pt2'], element['m_y_pt2'], element['m_z_pt2'])
            ele.set_vert(2, element['m_x_pt3'], element['m_y_pt3'], element['m_z_pt3'])
            ele.set_vert(3, element['m_x_pt4'], element['m_y_pt4'], element['m_z_pt4'])
            strikes.append(ele.strike())
            rakes.append(element['rake_rad'])
        
        self.look_azimuth = -sum(strikes)/len(strikes)
        
        average_rake = sum(rakes)/len(rakes)
        if average_rake >= math.pi/2.0:
            average_rake = math.pi - average_rake
        self.look_elevation = abs(average_rake)
    
    #---------------------------------------------------------------------------
    # Returns a PIL image of the masked displacement map using the current
    # values of the displacements. This map can then be combined into a still
    # or used as part of an animation.
    #---------------------------------------------------------------------------
    def create_displacement_image(self, fringes=True):
        #-----------------------------------------------------------------------
        # Set all of the plotting properties
        #-----------------------------------------------------------------------
    
        # properties that are fringes dependent
        if fringes:
            cmap            = self.dmc['cmap_f']
            water_color     = self.dmc['water_color_f']
            boundary_color  = self.dmc['boundary_color_f']
        else:
            cmap            = self.dmc['cmap']
            water_color     = self.dmc['water_color']
            boundary_color  = self.dmc['boundary_color']
            
        # properties that are not fringes dependent
        land_color      = cmap(0)
        plot_resolution = self.dmc['plot_resolution']
        
        #-----------------------------------------------------------------------
        # Set the map dimensions
        #-----------------------------------------------------------------------
        mw = self.lons_1d.size
        mh = self.lats_1d.size
        mwi = mw/plot_resolution
        mhi = mh/plot_resolution
        
        #-----------------------------------------------------------------------
        # Fig1 is the background land and ocean.
        #-----------------------------------------------------------------------
        fig1 = mplt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
        self.m1.ax = fig1.add_axes((0,0,1,1))
        self.m1.drawmapboundary(
            color=boundary_color,
            linewidth=0,
            fill_color=water_color
        )
        self.m1.fillcontinents(
            color=land_color,
            lake_color=water_color
        )
        
        #-----------------------------------------------------------------------
        # Fig2 is the deformations.
        #-----------------------------------------------------------------------
        fig2 = mplt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
        self.m2.ax = fig2.add_axes((0,0,1,1))
        
        if self.look_azimuth is None:
            self.look_azimuth = 0.0
        if self.look_elevation is None:
            self.look_elevation = 0.0
        
        dMags = -self.dX * math.sin(self.look_azimuth) * math.cos(self.look_elevation) - self.dY * math.cos(self.look_azimuth) * math.cos(self.look_elevation) + self.dZ * math.sin(self.look_elevation)
        
        
        #prepare the colors for the plot
        dMags_transformed = self.m2.transform_scalar(dMags, self.lons_1d, self.lats_1d, self.lons_1d.size, self.lats_1d.size)
        
        #print 'start colors'
        if fringes:
            dMags_colors = np.empty((dMags_transformed.shape[0],dMags_transformed.shape[1],4))
            r,g,b,a = cmap(0)
            dMags_colors[:,:,0].fill(r)
            dMags_colors[:,:,1].fill(g)
            dMags_colors[:,:,2].fill(b)
            dMags_colors[:,:,3].fill(a)
            non_zeros = dMags_transformed.nonzero()
            for n,i in enumerate(non_zeros[0]):
                j = non_zeros[1][n]
                r,g,b,a = cmap(math.modf(abs(dMags_transformed[i,j])/self.wavelength)[0])
                dMags_colors[i, j, 0] = r
                dMags_colors[i, j, 1] = g
                dMags_colors[i, j, 2] = b
                dMags_colors[i, j, 3] = a
            '''
            it = np.nditer(dMags_transformed, flags=['multi_index'])
            while not it.finished:
                r,g,b,a = cmap(math.modf(abs(dMags_transformed[it.multi_index])/self.wavelength)[0])
                dMags_colors[it.multi_index[0], it.multi_index[1], 0] = r
                dMags_colors[it.multi_index[0], it.multi_index[1], 1] = g
                dMags_colors[it.multi_index[0], it.multi_index[1], 2] = b
                dMags_colors[it.multi_index[0], it.multi_index[1], 3] = a
                it.iternext()
            #cPickle.dump(dMags_colors, open('local/test_colors.pkl','wb'))
            '''
            #dMags_colors = cPickle.load(open('local/test_colors.pkl','rb'))
            #print dMags_colors
            im = self.m2.imshow(dMags_colors, interpolation='spline36')
        else:
            dMags_colors = np.fabs(dMags_transformed)
            vmax = np.amax(dMags_colors)
            if vmax <= 1:
                mod_vmax = 1
            elif vmax > 1 and vmax <= 10:
                mod_vmax = 10
            elif vmax > 10 and vmax <= 100:
                mod_vmax = 100
            elif vmax > 100 and vmax <= 1000:
                mod_vmax = 1000
            elif vmax > 1000:
                mod_vmax = 1000
            im = self.m2.imshow(dMags_colors, cmap=cmap, norm=mcolor.LogNorm(vmin=1e-4, vmax=mod_vmax, clip=True))
        #print 'end colors'

        #-----------------------------------------------------------------------
        # Fig3 is the land/sea mask.
        #-----------------------------------------------------------------------
        fig3 = mplt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
        self.m3.ax = fig3.add_axes((0,0,1,1))
        self.m3.fillcontinents(color='#000000', lake_color='#ffffff')
        
        #-----------------------------------------------------------------------
        # Composite fig 1 - 3 together
        #-----------------------------------------------------------------------
        # FIGURE 1 draw the renderer
        fig1.canvas.draw()
        
        # FIGURE 1 Get the RGBA buffer from the figure
        w,h = fig1.canvas.get_width_height()
        buf = np.fromstring ( fig1.canvas.tostring_argb(), dtype=np.uint8 )
        buf.shape = ( w, h,4 )
     
        # FIGURE 1 canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
        buf = np.roll ( buf, 3, axis = 2 )
        im1 = Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
        
        # FIGURE 2 draw the renderer
        fig2.canvas.draw()
        
        # FIGURE 2 Get the RGBA buffer from the figure
        w,h = fig2.canvas.get_width_height()
        buf = np.fromstring ( fig2.canvas.tostring_argb(), dtype=np.uint8 )
        buf.shape = ( w, h,4 )
     
        # FIGURE 2 canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
        buf = np.roll ( buf, 3, axis = 2 )
        im2 = Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
        
        # FIGURE 3 draw the renderer
        fig3.canvas.draw()
        
        # FIGURE 3 Get the RGBA buffer from the figure
        w,h = fig3.canvas.get_width_height()
        buf = np.fromstring ( fig3.canvas.tostring_argb(), dtype=np.uint8 )
        buf.shape = ( w, h,4 )
     
        # FIGURE 3 canvas.tostring_argb give pixmap in ARGB mode. Roll the ALPHA channel to have it in RGBA mode
        buf = np.roll ( buf, 3, axis = 2 )
        im3 = Image.fromstring( "RGBA", ( w ,h ), buf.tostring( ) )
        
        mask = im3.convert('L')
        
        # The final composited image.
        return  Image.composite(im1, im2, mask)



#-------------------------------------------------------------------------------
# plots event displacements
#-------------------------------------------------------------------------------
def plot_event_displacements(sim_file, output_file, evnum, fringes=True, padding=0.01, cutoff=None):
    start_time = time.time()
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)
        
        min_lat = geometry.min_lat
        max_lat = geometry.max_lat
        min_lon = geometry.min_lon
        max_lon = geometry.max_lon
        base_lat = geometry.base_lat
        base_lon = geometry.base_lon

        event_data = events[evnum]
        event_element_slips = events.get_event_element_slips(evnum)
        ele_getter = itemgetter(*event_element_slips.keys())
        event_element_data = ele_getter(geometry)
        if len(event_element_slips) == 1:
            event_element_data = [event_element_data]
        fault_traces = geometry.get_fault_traces()
        event_sections = geometry.sections_with_elements(event_element_slips.keys())

    #print event_element_data
    print 'Done initilizing data - {} seconds'.format(time.time() - start_time)
    print '{} elements in event'.format(len(event_element_slips))
    print '{} sections in event'.format(len(event_sections))
    print event_sections
    
    start_time = time.time()
    dmp = VCDisplacementMapPlotter(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, fault_traces, padding=0.01)
    print 'Done initilizing plotter - {} seconds'.format(time.time() - start_time)

    start_time = time.time()
    dmp.calculate_displacements(event_element_data, event_element_slips, cutoff=cutoff)
    #dmp.dX, dmp.dY, dmp.dZ = cPickle.load(open('local/test_disp.pkl','rb'))
    print 'Done calculating displacements - {} seconds'.format(time.time() - start_time)
    cPickle.dump((dmp.dX, dmp.dY, dmp.dZ), open('local/test_disp.pkl','wb'))
    
    start_time = time.time()
    dmp.calculate_look_angles(event_element_data)
    print 'Done calculating look angles - {} seconds'.format(time.time() - start_time)
    
    start_time = time.time()
    map_image = dmp.create_displacement_image(fringes=fringes)
    print 'Done creating the map image - {} seconds'.format(time.time() - start_time)
    
    start_time = time.time()
    #---------------------------------------------------------------------------
    # Plot all of the geographic info on top of the displacement map image.
    #---------------------------------------------------------------------------
    
    # Grab all of the plot properties that we will need.
    # properties that are fringes dependent
    if fringes:
        cmap            = dmp.dmc['cmap_f']
        coastline_color = dmp.dmc['coastline_color_f']
        country_color   = dmp.dmc['country_color_f']
        state_color     = dmp.dmc['state_color_f']
        fault_color     = dmp.dmc['fault_color_f']
        map_tick_color  = dmp.dmc['map_tick_color_f']
        map_frame_color = dmp.dmc['map_frame_color_f']
        grid_color      = dmp.dmc['grid_color_f']
        cb_fontcolor    = dmp.dmc['cb_fontcolor_f']
    else:
        cmap            = dmp.dmc['cmap']
        coastline_color = dmp.dmc['coastline_color']
        country_color   = dmp.dmc['country_color']
        state_color     = dmp.dmc['state_color']
        fault_color     = dmp.dmc['fault_color']
        map_tick_color  = dmp.dmc['map_tick_color']
        map_frame_color = dmp.dmc['map_frame_color']
        grid_color      = dmp.dmc['grid_color']
        cb_fontcolor    = dmp.dmc['cb_fontcolor']
    
    # properties that are not fringes dependent
    boundary_width  = dmp.dmc['boundary_width']
    coastline_width = dmp.dmc['coastline_width']
    country_width   = dmp.dmc['country_width']
    state_width     = dmp.dmc['state_width']
    river_width     = dmp.dmc['river_width']
    fault_width     = dmp.dmc['fault_width']
    map_resolution  = dmp.dmc['map_resolution']
    plot_resolution = dmp.dmc['plot_resolution']
    map_frame_width = dmp.dmc['map_frame_width']
    map_fontsize    = dmp.dmc['map_fontsize']
    arrow_inset     = dmp.dmc['arrow_inset']
    arrow_fontsize  = dmp.dmc['arrow_fontsize']
    cb_fontsize     = dmp.dmc['cb_fontsize']
    cb_height       = dmp.dmc['cb_height']
    cb_margin_t     = dmp.dmc['cb_margin_t']
    grid_width      = dmp.dmc['grid_width']
    num_grid_lines  = dmp.dmc['num_grid_lines']
    font            = dmp.dmc['font']
    font_bold       = dmp.dmc['font_bold']
    map_resolution  = dmp.dmc['map_resolution']
    map_projection  = dmp.dmc['map_projection']
    
    # The sizing for the image is tricky. The aspect ratio of the plot is fixed,
    # so we cant set all of margins to whatever we want. We will set the anchor
    # to the top, left margin position. Then scale the image based on the
    # bottom/right margin, whichever is bigger.
    
    mw = dmp.lons_1d.size
    mh = dmp.lats_1d.size
    
    if mh > mw:
        ph = 768.0
        pw = mw + 70.0 + 40.0
    else:
        pw = 790.0
        ph = mh + 70.0 + 40.0

    width_frac = mw/pw
    height_frac = mh/ph
    left_frac = 70.0/pw
    bottom_frac = 70.0/ph

    pwi = pw/plot_resolution
    phi = ph/plot_resolution
    fig4 = mplt.figure(figsize=(pwi, phi), dpi=plot_resolution)

    #---------------------------------------------------------------------------
    # m4, fig4 is all of the boundary data.
    #---------------------------------------------------------------------------
    m4 = Basemap(
        llcrnrlon=dmp.padded_min_lon,
        llcrnrlat=dmp.padded_min_lat,
        urcrnrlon=dmp.padded_max_lon,
        urcrnrlat=dmp.padded_max_lat,
        lat_0=(dmp.padded_max_lat+dmp.padded_min_lat)/2.0,
        lon_0=(dmp.padded_max_lon+dmp.padded_min_lon)/2.0,
        resolution=map_resolution,
        projection=map_projection,
        suppress_ticks=True
    )
    m4.ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
    
    # draw coastlines, edge of map.
    m4.drawcoastlines(color=coastline_color, linewidth=coastline_width)
    
    # draw countries
    m4.drawcountries(linewidth=country_width, color=country_color)
    
    # draw states
    m4.drawstates(linewidth=state_width, color=state_color)
    
    # draw parallels.
    parallels = np.linspace(dmp.lats_1d.min(), dmp.lats_1d.max(), num_grid_lines+1)
    m4_parallels = m4.drawparallels(parallels, labels=[1,0,0,0], fontsize=map_fontsize, color=grid_color, fontproperties=font, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])
    
    # draw meridians
    meridians = np.linspace(dmp.lons_1d.min(), dmp.lons_1d.max(), num_grid_lines+1)
    m4_meridians = m4.drawmeridians(meridians, labels=[0,0,1,0], fontsize=map_fontsize, color=grid_color, fontproperties=font, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])
    
    # draw the azimuth look arrow
    az_width_frac    = 50.0/pw
    az_height_frac   = 50.0/ph
    az_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
    az_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac)/ph
    az_ax = fig4.add_axes((az_left_frac,az_bottom_frac,az_width_frac,az_height_frac))

    az_ax.set_xlim((0,1.0))
    az_ax.set_ylim((0,1.0))
    for item in az_ax.yaxis.get_ticklabels() + az_ax.xaxis.get_ticklabels() + az_ax.yaxis.get_ticklines() + az_ax.xaxis.get_ticklines():
        item.set_alpha(0)

    az_arrow_start_x    = 0.5 - (0.8/2.0)*math.sin(dmp.look_azimuth)
    az_arrow_start_y    = 0.5 - (0.8/2.0)*math.cos(dmp.look_azimuth)
    az_arrow_dx      = 0.8*math.sin(dmp.look_azimuth)
    az_arrow_dy      = 0.8*math.cos(dmp.look_azimuth)

    az_ax.arrow( az_arrow_start_x , az_arrow_start_y, az_arrow_dx, az_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='right', length_includes_head=True, lw=1.0, fc='k' )
    az_ax.add_line(mlines.Line2D((0.5,0.5), (0.5,0.8), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
    az_ax.add_patch(mpatches.Arc((0.5,0.5), 0.3, 0.3, theta1=90.0 - dmp.convert.rad2deg(dmp.look_azimuth), theta2=90.0, fc='none', lw=1.0, ls='dotted', ec='k'))
    az_ax.text(1.0, 1.0, 'az = {:0.1f}{}'.format(dmp.convert.rad2deg(dmp.look_azimuth),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')

    # draw the altitude look arrow
    al_width_frac    = 50.0/pw
    al_height_frac   = 50.0/ph
    al_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
    al_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac - ph*al_height_frac)/ph
    al_ax = fig4.add_axes((al_left_frac,al_bottom_frac,al_width_frac,al_height_frac))

    al_ax.set_xlim((0,1.0))
    al_ax.set_ylim((0,1.0))
    for item in al_ax.yaxis.get_ticklabels() + al_ax.xaxis.get_ticklabels() + al_ax.yaxis.get_ticklines() + al_ax.xaxis.get_ticklines():
        item.set_alpha(0)

    al_arrow_start_x    = 0.1 + 0.8*math.cos(dmp.look_elevation)
    al_arrow_start_y    = 0.1 + 0.8*math.sin(dmp.look_elevation)
    al_arrow_dx      = -0.8*math.cos(dmp.look_elevation)
    al_arrow_dy      = -0.8*math.sin(dmp.look_elevation)

    al_ax.arrow( al_arrow_start_x , al_arrow_start_y, al_arrow_dx, al_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='left', length_includes_head=True, lw=1.0, fc='k' )
    al_ax.add_line(mlines.Line2D((0.1,0.9), (0.1,0.1), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
    al_ax.add_patch(mpatches.Arc((0.1,0.1), 0.5, 0.5, theta1=0.0, theta2=dmp.convert.rad2deg(dmp.look_elevation), fc='none', lw=1.0, ls='dotted', ec='k'))
    al_ax.text(1.0, 1.0, 'al = {:0.1f}{}'.format(dmp.convert.rad2deg(dmp.look_elevation),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')
    
    # draw the box with the magnitude
    mag_width_frac    = 50.0/pw
    mag_height_frac   = 10.0/ph
    mag_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
    mag_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac  - ph*az_height_frac - ph*mag_height_frac)/ph
    mag_ax = fig4.add_axes((mag_left_frac,mag_bottom_frac,mag_width_frac,mag_height_frac))

    mag_ax.set_xlim((0,1.0))
    mag_ax.set_ylim((0,1.0))
    for item in mag_ax.yaxis.get_ticklabels() + mag_ax.xaxis.get_ticklabels() + mag_ax.yaxis.get_ticklines() + mag_ax.xaxis.get_ticklines():
        item.set_alpha(0)
    
    mag_ax.text(0.5, 0.5, 'm = {:0.3f}'.format(float(event_data['event_magnitude'])), fontproperties=font_bold, size=arrow_fontsize, ha='center', va='center')

    # add the displacement map image to the plot
    m4.imshow(map_image, origin='upper')
    
    # print faults on lon-lat plot
    for sid, sec_trace in dmp.fault_traces.iteritems():
        trace_Xs, trace_Ys = m4(sec_trace[1], sec_trace[0])
        
        if sid in event_sections:
            linewidth = fault_width + 3
        else:
            linewidth = fault_width

        m4.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=linewidth, solid_capstyle='round', solid_joinstyle='round')

    #plot the cb
    cb_left_frac = 70.0/pw
    cb_bottom_frac = (70.0 - cb_height - cb_margin_t)/ph
    cb_width_frac = width_frac
    cb_height_frac = cb_height/ph

    cb_ax = fig4.add_axes([cb_left_frac, cb_bottom_frac, cb_width_frac, cb_height_frac])
    if fringes:
        norm = mcolor.Normalize(vmin=0, vmax=dmp.wavelength)
        cb = mcolorbar.ColorbarBase(cb_ax, cmap=cmap,
                   norm=norm,
                   orientation='horizontal')
        cb_ax.set_title('Displacement [m]', fontproperties=font, color=cb_fontcolor, size=cb_fontsize, va='top', ha='left', position=(0,-1.5) )
    else:
        norm = mcolor.LogNorm(vmin=1e-4, vmax=mod_vmax, clip=True)
        cb = mcolorbar.ColorbarBase(cb_ax, cmap=cmap,
                   norm=norm,
                   orientation='horizontal')
        cb_ax.set_title('Total displacement [m]', fontproperties=font, color=cb_fontcolor, size=cb_fontsize, va='top', ha='left', position=(0,-1.5) )

    for label in cb_ax.xaxis.get_ticklabels():
        label.set_fontproperties(font)
        label.set_fontsize(cb_fontsize)
        label.set_color(cb_fontcolor)
    for line in cb_ax.xaxis.get_ticklines():
        line.set_alpha(0)

    # save the figure
    fig4.savefig(output_file, format='png', dpi=plot_resolution)

    print 'Done plotting the image - {} seconds'.format(time.time() - start_time)

#-------------------------------------------------------------------------------
# plots recurrence intervals
#-------------------------------------------------------------------------------
def plot_recurrence_intervals(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None):
    #---------------------------------------------------------------------------
    # Plot setup
    #---------------------------------------------------------------------------
    
    num_cols = 5.0
    
    # dimensions
    simw = 270.0
    simh = 270.0
    stm = 40.0
    sbm = 40.0
    slm = 50.0
    srm = 10.0
    res = 72.0
    
    # fonts
    ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
    framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=10)
    legendfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
    titlefont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=12)
    subtitlefont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=8)
    
    #---------------------------------------------------------------------------
    # Instantiate the VCSimData class using the with statement. Then instantiate
    # VCEvents class from within the with block. This ensures that the sim data
    # file is closed when the with block ends.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc events class passing in an instance of the
        # VCSimData class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)
        
        # get the data
        section_info = geometry.get_section_info(section_filter=section_filter)
        
        #-----------------------------------------------------------------------
        # start the plot
        #-----------------------------------------------------------------------
        bins = np.linspace(0,250,50)
        # calculate the final dimensions and create the figure and axis
        spw = simw - slm - srm
        sph = simh - stm - sbm
        num_rows = math.ceil(float(len(section_info))/num_cols)
        imw = math.ceil(simw * num_cols)
        imh = math.ceil(simh * num_rows)
        imwi = imw/res
        imhi = imh/res
        fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
        
        #-----------------------------------------------------------------------
        # Calculate the recurrence intervals and plot.
        #-----------------------------------------------------------------------
        curr_row = -1.0
        for num, secid in enumerate(sorted(section_info.keys())):
            curr_col = num%num_cols
            if curr_col == 0.0:
                curr_row += 1.0
            #print curr_row, curr_col
            the_ax = fig.add_axes(((slm + curr_col * simw)/imw, (sbm + (num_rows - curr_row - 1) * simh)/imh, spw/imw, sph/imh))
            section_events = events.get_event_data_from_evids(
                                        geometry.events_on_section(secid),
                                        ['event_magnitude', 'event_year'],
                                        event_range=event_range,
                                        magnitude_filter='>=6.5'
                                    )
            intervals = [
                x - section_events['event_year'][n-1]
                for n,x in enumerate(section_events['event_year'])
                if n != 0]
                
            intervals7 = [
                x - section_events['event_year'][n-1]
                for n,x in enumerate(section_events['event_year'])
                if n != 0 and section_events['event_magnitude'][n] >= 7.0]
            
            
            hist, bins = np.histogram(intervals, bins=bins, density=True)
            hist7, bins7 = np.histogram(intervals7, bins=bins, density=True)
            mean = np.mean(intervals)
            std = np.std(intervals)
            mean7 = np.mean(intervals7)
            std7 = np.std(intervals7)
            
            the_ax.step(bins[0:-1], hist, where='post', label='m>6.5')
            the_ax.step(bins7[0:-1], hist7, where='post', label='m>7')
            
            for label in the_ax.xaxis.get_ticklabels()+the_ax.yaxis.get_ticklabels():
                label.set_fontproperties(ticklabelfont)
                
            the_ax.set_ylabel('Prob. Density', fontproperties=framelabelfont)
            the_ax.set_xlabel('Recurrence Time [yr]', fontproperties=framelabelfont)
            
            the_ax.autoscale_view(tight=True)
            
            the_ax.set_title('{} {}'.format(secid,section_info[secid]['name']), position=(0.0,1.04), ha='left', fontproperties=titlefont)
            the_ax.text(0.0, 1.01, 'm>6.5: mean {mean:0.1f} std {std:0.1f}, m>7.0 mean {mean7:0.1f} std {std7:0.1f}'.format(mean=mean, std=std, mean7=mean7, std7=std7), va='bottom', ha='left', transform=the_ax.transAxes, fontproperties=subtitlefont)
    
            the_ax.legend(prop=legendfont)

    # Get the plot format and save the file
    plot_format = output_file.split('.')[-1]
    if plot_format != 'png' and plot_format != 'pdf':
        raise vcexceptions.PlotFormatNotSupported(plot_format)
    else:
        fig.savefig(output_file, format=plot_format, dpi=res)

#-------------------------------------------------------------------------------
# plots an event graph
#-------------------------------------------------------------------------------
def plot_graph(graph_file, output_file, degree_cut=None, label_degree_cut=0.25, self_loops=True):
    G = cPickle.load(open(graph_file, 'rb'))
    
    #print(nx.clustering(nx.Graph(G), weight='weight'))
    
    # the color map for the plot
    cmap = mplt.get_cmap('GnBu_r')
    
    if degree_cut is not None:
        print 'Original Graph'
        print nx.info(G)
        degrees = G.degree(weight='weight')
        max_degree = float(max(degrees.values()))
        min_degree = float(min(degrees.values()))
        degree_cut_num = min_degree + (max_degree-min_degree)*degree_cut
        print 'max degree: {}'.format(max_degree)
        print 'min degree: {}'.format(min_degree)
        print 'degree cut: {}'.format(degree_cut_num)
        print
        print 'Cut Graph'
        sub_nodes = [n for n, d in G.degree(weight='weight').iteritems() if d > degree_cut_num]
        Gsub = G.subgraph(sub_nodes)
    else:
        Gsub = G

    if not self_loops:
        print 'Removing Self Loops'
        self_loop_edges = Gsub.selfloop_edges()
        Gsub.remove_edges_from(self_loop_edges)
    
    print nx.info(Gsub)
    
    degrees = Gsub.degree(weight='weight')
    max_degree = float(max(degrees.values()))
    min_degree = float(min(degrees.values()))
    print 'max degree: {}'.format(max_degree)
    print 'min degree: {}'.format(min_degree)
    node_min = 0.01
    node_max = 0.2
    node_line_min = 0.1
    node_line_max = 2.0
    min_label_degree = min_degree + (max_degree-min_degree)*label_degree_cut
    min_font_size = 0.5
    max_font_size = 6.0

    widths = {}
    heights = {}
    labels = {}
    styles = {}
    colors = {}
    node_line_widths = {}
    font_sizes = {}
    #print max_degree, min_degree
    for n in nx.nodes_iter(Gsub):
        degree = float(Gsub.degree(n, weight='weight'))
        r,g,b,a = cmap(vcutils.linear_interp(degree, min_degree, max_degree, 0.0, 1.0))
        dim = vcutils.linear_interp(degree, min_degree, max_degree, node_min, node_max)
        widths[n] = dim
        heights[n] = dim
        if degree > min_label_degree:
            labels[n] = n
            font_sizes[n] = vcutils.linear_interp(degree, min_degree, max_degree, min_font_size, max_font_size)
        else:
            labels[n] = ''
        styles[n] = 'filled'
        colors[n] = '#{r:02x}{g:02x}{b:02x}'.format(r=int(r*255.0), g=int(g*255.0), b=int(b*255.0))
        node_line_widths[n] = vcutils.linear_interp(degree, min_degree, max_degree, node_line_min, node_line_max)

    nx.set_node_attributes(Gsub,'width',widths)
    nx.set_node_attributes(Gsub,'height',heights)
    nx.set_node_attributes(Gsub,'label',labels)
    nx.set_node_attributes(Gsub,'style',styles)
    nx.set_node_attributes(Gsub,'fillcolor',colors)
    nx.set_node_attributes(Gsub,'penwidth',node_line_widths)
    nx.set_node_attributes(Gsub,'fontsize',font_sizes)
    #print G.edges(data=True)
    
    weights = [ float(edata['weight']) for u,v,edata in Gsub.edges(data=True) ]
    
    max_weight = float(max(weights))
    min_weight = float(min(weights))
    line_min = 0.1
    line_max = 5.0
    
    edge_widths = {}
    arrow_sizes = {}
    edge_colors = {}
    for e in nx.edges_iter(Gsub):
        width = vcutils.linear_interp(float(Gsub[e[0]][e[1]]['weight']), min_weight, max_weight, line_min, line_max)
        alpha = vcutils.linear_interp(float(Gsub[e[0]][e[1]]['weight']), min_weight, max_weight, 10.0, 255.0)
        edge_widths[e] = width
        arrow_sizes[e] = 0.1
        edge_colors[e] = '#000000{:x}'.format(int(alpha))
    
    nx.set_edge_attributes(Gsub, 'penwidth', edge_widths)
    nx.set_edge_attributes(Gsub, 'arrowsize', arrow_sizes)
    nx.set_edge_attributes(Gsub, 'color', edge_colors)
    #cmap = mplt.get_cmap('gray')
    
    #norm = mcolor.Normalize(vmin=min(edge_weights), vmax=max(edge_weights))
    
    A=nx.to_agraph(Gsub)        # convert to a graphviz graph
    #A.layout()            # neato layout
    A.draw(output_file, prog='sfdp', args='-Gsize="40!" -Goverlap=prism -Grepulsiveforce=1.0 -GsmoothType="graph_dist" -Goutputorder="edgesfirst" -Nfixedsize="true" -Nfontname="Helvetica"')

#-------------------------------------------------------------------------------
# space-time plot
#-------------------------------------------------------------------------------
def space_time_plot(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None):
    #---------------------------------------------------------------------------
    # Instantiate the VCSimData class using the with statement. Then instantiate
    # VCEvents class from within the with block. This ensures that the sim data
    # file is closed when the with block ends.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)

        # get the data
        event_data = events.get_event_data(['event_number','event_year','event_magnitude','event_elements', 'event_range_duration'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        section_info = geometry.get_section_info(section_filter=section_filter)
        
        # store a sorted list of section ids
        section_ids = sorted(section_info.keys())
        
        # the section offsets determine the starting x position of each section
        section_offsets = {}
        for i, sid in enumerate(section_ids):
            section_offsets[sid] = sum([section_info[k]['blocks_along_strike'] for k in sorted(section_info.keys())[0:i]])
        
        # calculate various properties of the data set that we will need to
        # set up the plot
        min_depth = min([section_info[k]['blocks_along_dip'] for k in section_info.keys()])
        x_data_size = sum([section_info[k]['blocks_along_strike'] for k in section_info.keys()])
        max_label_len = max([len(section_info[k]['name']) for k in section_info.keys()])
        start_year = event_data['event_year'][0]
        
        # Storing all of the plot parameters here for clarity
        stp_params = {
            'output_file':output_file,
            'x_axis_data_size':x_data_size,
            'y_axis_data_size':event_data['event_range_duration'],
            'max_depth':min_depth,
            'min_mag':min(event_data['event_magnitude']),
            'max_mag':max(event_data['event_magnitude']),
            'start_year':start_year,
            'max_label_len':max_label_len,
            'geometry':geometry,
            'section_offsets':section_offsets
        }
        
        # instantiate the spacetimeplot class
        stp = vcutils.VCSpaceTimePlot(
            stp_params['output_file'],
            stp_params['x_axis_data_size'],
            stp_params['y_axis_data_size'],
            stp_params['max_depth'],
            stp_params['min_mag'],
            stp_params['max_mag'],
            stp_params['start_year'],
            stp_params['max_label_len']
        )
        
        mp = False
        #-----------------------------------------------------------------------
        # The multiprocessing stuff below is not functional. The variable "mp"
        # above should always be set to False.
        #-----------------------------------------------------------------------
        # TODO: Figure out a way to plot in parallel.
        if mp:
            num_processes = multiprocessing.cpu_count()
        
            # break the work up
            seg = int(round(float(len(event_data['event_magnitude']))/float(num_processes)))
            work_queue = multiprocessing.Queue()
            for i in range(num_processes):
                if i == num_processes - 1:
                    end_index = len(event_data['event_magnitude'])
                else:
                    end_index = seg*int(i + 1)
                work_queue.put({
                    'event_magnitude':event_data['event_magnitude'][int(i) * seg:end_index],
                    'event_elements':event_data['event_elements'][int(i) * seg:end_index],
                    'event_number':event_data['event_number'][int(i) * seg:end_index],
                    'event_year':event_data['event_year'][int(i) * seg:end_index]
                })

            # create a queue to pass to workers to store the results
            result_queue = multiprocessing.Queue()

            # spawn workers
            for i in range(num_processes):
                worker = vcutils.SpaceTimePlotter(stp_params, work_queue, result_queue)
                worker.start()
            
            # collect the results off the queue
            for i in range(num_processes):
                stp.event_lines += result_queue.get().event_lines
        else:
            # For each event in the found event set, look at the involved
            # elements, and add them to the event line array. Since the event
            # line shows only elements on the strike, elements at depths are
            # projected up to the strike: for every element along the dip the
            # strike value is incremented up to the smallest value of depth in
            # the model.
            for i, enum in enumerate(event_data['event_number']):
                event_line = np.zeros(x_data_size)
                for bid in event_data['event_elements'][i]:
                    sid = geometry[bid]['section_id']
                    try:
                        b_index = section_offsets[sid] + geometry[bid]['das_id']
                        if event_line[b_index] < min_depth:
                            event_line[b_index] += 1
                    except KeyError:
                        pass
                stp.add_event(
                    enum,
                    event_data['event_year'][i],
                    event_data['event_magnitude'][i],
                    event_line
                )

        # Add section labels
        stp.add_section_labels(section_offsets, section_info)

        # Add the title
        stp.add_title('Events from {}'.format(sim_data.filename))

        # Plot the thing
        stp.plot()
    
#-------------------------------------------------------------------------------
# magnitude rupture area plot
#-------------------------------------------------------------------------------
def magnitude_rupture_area(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None):
    #---------------------------------------------------------------------------
    # Instantiate the VCSimData class using the with statement. Then instantiate
    # VCEvents class from within the with block. This ensures that the sim data
    # file is closed when the with block ends.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc events class passing in an instance of the
        # VCSimData class
        events = VCEvents(sim_data)
        
        # get the data
        event_data = events.get_event_data(['event_magnitude', 'event_area'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
    
    #---------------------------------------------------------------------------
    # Prepare the plot and do it.
    #---------------------------------------------------------------------------
    # TODO: Move this to another function
    
    # All of the data is in mks units. We need kilometers for this plot.
    event_area_kmsq = [vcutils.Converter().msq_kmsq(x) for x in event_data['event_area']]
    
    # get the binned averages of the data
    x_ave, y_ave = vcutils.calculate_averages(event_area_kmsq, event_data['event_magnitude'])
    
    # get the plot label which will depend on the filters
    plot_label = vcutils.get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)
    
    x_WC = np.linspace(2.2,5184)
    y_WC = 4.07 + 0.98 * np.log10(x_WC)
    y_error_plus_WC = 4.07+0.06 + (0.98+0.03) * np.log10(x_WC)
    y_error_minus_WC = 4.07-0.06 + (0.98-0.03) * np.log10(x_WC)
    y_error_WC = [np.subtract(y_WC, y_error_minus_WC), np.subtract(y_error_plus_WC, y_WC)]

    # do the standard plot
    vcutils.standard_plot(output_file, event_area_kmsq, event_data['event_magnitude'],
        axis_format='semilogx',
        add_lines=[
            {'label':'binned average', 'x':x_ave, 'y':y_ave},
            {'label':'WC', 'x':x_WC, 'y':y_WC, 'ls':'--', 'c':'red'}
        ],
        axis_labels = {'x':r'log(Rupture Area [km$^\mathsf{2}$])', 'y':'Magnitude'},
        plot_label='Magnitude-Rupture Area{}'.format(plot_label)
    )
        
#-------------------------------------------------------------------------------
# magnitude average slip plot
#-------------------------------------------------------------------------------
def magnitude_average_slip(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None):
    #---------------------------------------------------------------------------
    # Instantiate the VCSimData class using the with statement. Then instantiate
    # VCEvents class from within the with block. This ensures that the sim data
    # file is closed when the with block ends.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc events class passing in an instance of the
        # VCSimData class
        events = VCEvents(sim_data)
        
        # get the data
        event_data = events.get_event_data(['event_magnitude', 'event_average_slip'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
    
    #---------------------------------------------------------------------------
    # Prepare the plot and do it.
    #---------------------------------------------------------------------------
    # TODO: Move this to another function

    # get the binned averages of the data
    x_ave, y_ave = vcutils.calculate_averages(event_data['event_average_slip'], event_data['event_magnitude'])
    
    # get the plot label which will depend on the filters
    plot_label = vcutils.get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)

    x_WC = np.linspace(0.05,8, num=10)
    y_WC = 6.93 + 0.82 * np.log10(x_WC)
    y_error_plus_WC = 6.93+0.05 + (0.82+0.1) * np.log10(x_WC)
    y_error_minus_WC = 6.93-0.05 + (0.82-0.1) * np.log10(x_WC)
    y_error_WC = [np.subtract(y_WC, y_error_minus_WC), np.subtract(y_error_plus_WC, y_WC)]
    
    # do the standard plot
    vcutils.standard_plot(output_file, event_data['event_average_slip'], event_data['event_magnitude'],
        axis_format='semilogx',
        add_lines=[
            {'label':'binned average', 'x':x_ave, 'y':y_ave},
            {'label':'WC', 'x':x_WC, 'y':y_WC, 'ls':'--', 'c':'red'}
        ],
        axis_labels = {'y':'Magnitude', 'x':'log(Average Slip [m])'},
        plot_label='Magnitude-Average Slip{}'.format(plot_label)
    )

#-------------------------------------------------------------------------------
# average slip surface rupture length plot
#-------------------------------------------------------------------------------
def average_slip_surface_rupture_length(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None):
    #---------------------------------------------------------------------------
    # Instantiate the VCSimData class using the with statement. Then instantiate
    # VCEvents class from within the with block. This ensures that the sim data
    # file is closed when the with block ends.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc events class passing in an instance of the
        # VCSimData class
        events = VCEvents(sim_data)
        
        # get the data
        event_data = events.get_event_data(['event_surface_rupture_length', 'event_average_slip'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
    
    #---------------------------------------------------------------------------
    # Prepare the plot and do it.
    #---------------------------------------------------------------------------
    # TODO: Move this to another function
    
    # All of the data is in mks units. We need kilometers for this plot.
    event_surface_rupture_length_km = [vcutils.Converter().m_km(x) for x in event_data['event_surface_rupture_length']]
    
    # get the binned averages of the data
    x_ave, y_ave = vcutils.calculate_averages(event_surface_rupture_length_km, event_data['event_average_slip'])
    
    # get the plot label which will depend on the filters
    plot_label = vcutils.get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)

    x_WC = np.linspace(3.8,432, num=10)
    y_WC = 10.0**(-1.43 + 0.88 * np.log10(x_WC))
    y_error_plus_WC = 10.0**(-1.43+0.18 + (0.88+0.11) * np.log10(x_WC))
    y_error_minus_WC = 10.0**(-1.43-0.18 + (0.88-0.11) * np.log10(x_WC))
    y_error_WC = [np.subtract(y_WC, y_error_minus_WC), np.subtract(y_error_plus_WC, y_WC)]
    
    # do the standard plot
    vcutils.standard_plot(output_file, event_surface_rupture_length_km, event_data['event_average_slip'],
        axis_format='loglog',
        add_lines=[
            {'label':'binned average', 'x':x_ave, 'y':y_ave},
            {'label':'WC', 'x':x_WC, 'y':y_WC, 'ls':'--', 'c':'red'}
        ],
        axis_labels = {'y':'log(Average Slip [m])', 'x':'log(Surface Rupture Length [km])'},
        plot_label='Average Slip-Surface Rupture Length{}'.format(plot_label)
    )

#-------------------------------------------------------------------------------
# frequency magnitude plot
#-------------------------------------------------------------------------------
def frequency_magnitude(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None):
    #---------------------------------------------------------------------------
    # Instantiate the VCSimData class using the with statement. Then instantiate
    # VCEvents class from within the with block. This ensures that the sim data
    # file is closed when the with block ends.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc events class passing in an instance of the
        # VCSimData class
        events = VCEvents(sim_data)
        
        # get the data
        event_data = events.get_event_data(['event_magnitude', 'event_range_duration'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
    
    #---------------------------------------------------------------------------
    # Prepare the plot and do it.
    #---------------------------------------------------------------------------
    # TODO: Move this to another function
    
    # initilize a dict to store the event counts and get the total number
    # of events.
    cum_freq = {}
    total_events = len(event_data['event_magnitude'])
    
    # count the number of events bigger than each magnitude
    for num, magnitude in enumerate(sorted(event_data['event_magnitude'])):
        cum_freq[magnitude] = total_events - (num + 1)
    
    # dump the counts into x and y arrays for plotting. also, divide the count
    # by the number of years so we get count per year.
    x = []
    y = []
    for magnitude in sorted(cum_freq.iterkeys()):
        x.append(magnitude)
        y.append(float(cum_freq[magnitude])/event_data['event_range_duration'])

    # create the line for b = 1
    x_b1 = np.linspace(min(x),max(x),10)
    y_b1 = 10**(math.log(y[0],10)+x[0]-x_b1)

    # get the plot label which will depend on the filters
    plot_label = vcutils.get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)
    
    # for the UCERF2 error bars
    x_UCERF = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
    y_UCERF = [4.73, 2.15, 0.71, 0.24, 0.074, 0.020]
    y_error_UCERF = [[1.2, 0.37, 0.22, 0.09, 0.04, 0.016],[1.50, 0.43, 0.28, 0.11, 0.06, 0.035]]
    
    # do the standard plot
    vcutils.standard_plot(output_file, x, y,
        axis_format='semilogy',
        add_lines=[{'label':'b=1', 'x':x_b1, 'y':y_b1}, {'label':'UCERF2', 'x':x_UCERF, 'y':y_UCERF, 'ls':'--', 'c':'red', 'y_error':y_error_UCERF}],
        axis_labels = {'y':'log(# events per year)', 'x':'Magnitude'},
        plot_label='Frequency-Magnitude{}'.format(plot_label),
        connect_points=True,
        legend_loc='upper right'
    )