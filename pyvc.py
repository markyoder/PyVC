#!/usr/bin/env python

import site
import sys

site.addsitedir('')

#for i in sys.path:
#    print i

from optparse import OptionParser
import os

#VC specific modules
#sys.path.append('%s/classes/'%sys.path[0])
#sys.path.append('/Users/sachs/Documents/VirtualCalifornia/trunk/new_vcal/eqsim/')
#import Mesh
import VC
import Converter
import MyVector

import types
import cPickle

#from guppy import hpy

def get_refcounts():
    d = {}
    sys.modules
    # collect all classes
    for m in sys.modules.values():
        for sym in dir(m):
            o = getattr (m, sym)
            if type(o) is types.ClassType:
                d[o] = sys.getrefcount (o)
    # sort by refcount
    pairs = map (lambda x: (x[1],x[0]), d.items())
    pairs.sort()
    pairs.reverse()
    return pairs

def print_top_10():
    for n, c in get_refcounts()[:10]:
        print '%10d %s' % (n, c.__name__)  

def main(argv=None):
    if argv is None:
        argv = sys.argv
    
    parser = OptionParser()
    parser.add_option("-o", "--options-file", "--of",
                dest="options_file", default=None,
                help="THIS IS CURRENTLY DISABLED! The location of a file containing options for the script. If this option is set any other options sent in on the command line will be ignored.", 
                metavar="FILE"
                )
    parser.add_option("-d", "--data-file", "--df",
                dest="data_file", default=None,
                help="THIS CURRENTLY WORKS ONLY WITH HDF5 FILES! The location of a Virtual California hdf5 output file, a EqSim Geometry file, or an EqSim Events file. The type of file depends on the actions chosen. Some actions will only work with certian types of files.", 
                metavar="FILE"
                )
    parser.add_option("-t", "--data-file-type", "--dft",
                dest="data_file_type", default='h5', choices=('h5', 'eqsimg', 'eqsime'),
                help="THIS CURRENTLY WORKS ONLY WITH HDF5 FILES SO IT DOES NOT NEED TO BE SET! The type of data file we are parsing. The allowed values are h5, eqsimg and, eqsime. If it is an EqSim Geometry/Events file the conjugate Events/Geometry file will be needed. The directory where the data file is will be searched for this file. If no matching file can be found the user will be prompted to enter the path to the file."
                )
    parser.add_option("-f", "--fault-trace", "--ft",
                dest="fault_trace", default=None,
                help="The location of a fault trace file from which to create new geometry. If used in combination with a data file it will add the new geometry to the loaded geometry. The system name will be based on the name of the fault trace file.", 
                metavar="FILE"
                )
    parser.add_option("-s", "--sys-name",
                dest="sys_name", default=None,
                help="The system name. If not supplied the system name will be set based on the name of the data file or the fault trace file.", 
                metavar="SYS_NAME"
                )
    parser.add_option("-b", "--debug-mode",
                dest="debug_mode", action="store_true", default=False,
                help="Activate debug mode. Stops the system from caching certian stuff.", 
                )
    
    #These are filters for actions
    parser.add_option("--section-filter", "--sf",
                dest="section_filter", default=None,
                help="A comma seperated string of section IDs to include in any analysis performed.", 
                metavar="SEC_ID,...,SEC_ID"
                )
    parser.add_option("--start-year", "--sy",
                dest="start_year", type="int", default=None,
                help="The YEAR at which to start loading events. If not set this defaults to the first year in the data file. If used in conjunction with --start-event the value of --start-year will be used.", 
                metavar="YEAR"
                )
    parser.add_option("--end-year", "--ey",
                dest="end_year", type="int", default=None,
                help="The YEAR at which to stop loading events. If not set this defaults to the last year in the data file. If used in conjunction with --end-event the value of --end-year will be used.", 
                metavar="YEAR"
                )
    parser.add_option("--start-event", "--se",
                dest="start_event", type="int", default=0,
                help="The EVENT_ID at which to start loading events. If not set this defaults to the first event in the data file.", 
                metavar="EVENT_ID"
                )
    parser.add_option("--end-event", "--ee",
                dest="end_event", type="int", default=None,
                help="The EVENT_ID at which to stop loading events. If not set this defaults to the last event in the data file.", 
                metavar="EVENT_ID"
                )
    parser.add_option("--area-filter", "--af",
                dest="area_filter", default=None,
                help="A comma seperated string lat_NE,lon_NE,lat_SW,lon_SW or full for the full system", 
                metavar="SEC_ID,...,SEC_ID"
                )
    parser.add_option("--movie-only", "--mo",
                dest="movie_only", action="store_true",  default=False,
                help="If the frames of the displacement animation have already been created, setting this will just make the movie."
                )
    parser.add_option("--fringes",
                dest="fringes", action="store_true",  default=False,
                help="If set then the displacement map will plot interference fringes."
                )
    parser.add_option("--look-angles", "--la",
                dest="look_angles", default=None,
                help="A comma seperated string: look_azimuth [deg], look_elevation [deg]"
                )
    
    #These are the actions
    parser.add_option("--plot-faults", "--pf",
                dest="plot_faults", action="store_true",  default=False,
                help="Plot faults to a KML file."
                )
    parser.add_option("--single-event-animation", "--sea",
                dest="single_event_animation", type="int",  default=None,
                help="THIS IS CURRENTLY DISABLED! Single event animation of event EVENT_NUM.",
                metavar="EVENT_NUM"
                )
    parser.add_option("--export-eqsim-geometry", "--eeg",
                dest="export_eqsim_geometry", action="store_true",  default=False,
                help="Export the EQSim geometry files."
                )
    parser.add_option("--event-animation", "--ea",
                dest="event_animation", action="store_true",  default=False,
                help="Create an event animation in KML format."
                )
    parser.add_option("--rupture-map", "--rm",
                dest="rupture_map", action="store_true",  default=False,
                help="Create a rupture map."
                )
    parser.add_option("--frequency-magnitude", "--fm",
                dest="frequency_magnitude", choices=('both', 'text', 'plot'),  default=False,
                help="Create a frequency magnitude plot. OUTPUT_TYPE can be text, plot or both. In the text file that is output the columns are as follows: magnitude, number of greater than magnitude per year.",
                metavar="OUTPUT_TYPE"
                )
    parser.add_option("--magnitude-rupture-area", "--mra",
                dest="magnitude_rupture_area", choices=('both', 'text', 'plot'),  default=False,
                help="Create a magnitude rupture area plot. OUTPUT_TYPE can be text, plot or both. In the text file that is output the columns are as follows: rupture area [km^2], magnitude, binned average rupture area [km^2], binned average magnitude.",
                metavar="OUTPUT_TYPE"
                )
    parser.add_option("--magnitude-average-slip", "--mas",
                dest="magnitude_average_slip", choices=('both', 'text', 'plot'),  default=False,
                help="Create a magnitude average slip plot. OUTPUT_TYPE can be text, plot or both. In the text file that is output the columns are as follows: average slip [m], magnitude, binned average slip [m], binned magnitude.",
                metavar="OUTPUT_TYPE"
                )
    parser.add_option("--average-slip-rupture-length", "--asrl",
                dest="average_slip_rupture_length", choices=('both', 'text', 'plot'),  default=False,
                help="Create an average slip rupture length plot. OUTPUT_TYPE can be text, plot or both. In the text file that is output the columns are as follows: surface rupture length [km], average slip [m], binned surface rupture length [km], binned average slip [m].",
                metavar="OUTPUT_TYPE"
                )
    parser.add_option("--export-eqsim-events", "--eee",
                dest="export_eqsim_events", action="store_true",  default=False,
                help="Export eqsim events.",
                )
    parser.add_option("--eqsim-events-slip-map-min-mag", "--eesmmm",
                dest="eqsim_events_min_slip_map_mag", default=5.5,
                help="The minimum magnitude for the slip map in the EqSim event export. Defaults to 5.5.",
                metavar="MIN_SLIP_MAP_MAGNITUDE"
                )
    parser.add_option("--plot-format","--pltf",
                dest="plot_format", choices=('pdf', 'png'),  default='pdf',
                help="The format of the plot. Can be either pdf or png. Defaults to pdf.",
                metavar="PLOT_FORMAT"
                )
    parser.add_option("--element-CFF", "--cff",
                dest="element_CFF", choices=('both', 'text', 'plot'),  default=False,
                help="THIS IS EXPIRIMENTAL! Plot element CFF vs year. Only use this with small systems or use a section filter. OUTPUT_TYPE can be text, plot or both.",
                metavar="OUTPUT_TYPE"
                )
    parser.add_option("--element-CFF-map", "--cffmap",
                dest="element_CFF_map", action="store_true",  default=False,
                help="THIS IS EXPIRIMENTAL! Plot the map (in the dynamical systems sense) for each element.  Only use this with small systems or use a section filter.",
                )
    parser.add_option("--element-CFF-block-entropy", "--cffbe",
                dest="element_CFF_block_entropy", type="int",  default=None,
                help="THIS IS EXPIRIMENTAL! Plots the block entropy as a function of word length for each element. Set WORD_LENGTH to an integer to plot to a specific word length.  Only use this with small systems or use a section filter.",
                metavar="WORD_LENGTH"
                )
    parser.add_option("--element-CFF-delta-block-entropy", "--cffdbe",
                dest="element_CFF_delta_block_entropy", type="int",  default=None,
                help="THIS IS EXPIRIMENTAL! Plots the block entropy difference as a function of word length for each element. Set WORD_LENGTH to an integer to plot to a specific word length.  Only use this with small systems or use a section filter.",
                metavar="WORD_LENGTH"
                )
    parser.add_option("--element-CFF-E-Hmu", "--cffehmu",
                dest="element_CFF_E_Hmu", type="int",  default=None,
                help="THIS IS EXPIRIMENTAL! Plots the excess entropy and the entropy rate as a function of two dimensional position in the fault. Set WORD_LENGTH to an integer to plot to a specific word length.  Only use this with small systems or use a section filter.",
                metavar="WORD_LENGTH"
                )
    parser.add_option("--event-displacement-map", "--edm",
                dest="event_displacement_map", type="int",  default=None,
                help="Doc.",
                metavar="WORD_LENGTH"
                )
    parser.add_option("--event-displacement-map-all", "--edmall",
                dest="event_displacement_map_all", action="store_true",  default=False,
                help="Doc.",
                metavar="WORD_LENGTH"
                )
    parser.add_option("--event-displacement-map-animation", "--edma",
                dest="event_displacement_map_animation", action="store_true",  default=False,
                help="Doc.",
                metavar="WORD_LENGTH"
                )
    parser.add_option("--print-event-info", "--pei",
                dest="print_event_info", action="store_true",  default=False,
                help="Print event info.",
                )
    parser.add_option("--event-rupture-map", "--erm",
                dest="event_rupture_map", type="int",  default=None,
                help="Doc.",
                metavar="WORD_LENGTH"
                )
    parser.add_option("--event-fault-map", "--efm",
                dest="event_fault_map", type="int",  default=None,
                help="Doc.",
                metavar="WORD_LENGTH"
                )
    parser.add_option("--event-time-series", "--ets",
                dest="event_time_series", action="store_true",  default=False,
                help="Plot a time series of all the events in the simulation.",
                )
    parser.add_option("--event-time-series-marked", "--etsm",
                dest="event_time_series_marked", type="int",  default=None,
                help="Doc.",
                metavar="WORD_LENGTH"
                )
    parser.add_option("--plot-fault-map", "--pfm",
                dest="plot_fault_map", action="store_true",  default=False,
                help="Plot fault map to a png file. If a section filter is defined those sections will be marked on the map."
                )
    parser.add_option("--slip-in-slip-out", "--siso",
                dest="slip_in_slip_out", action="store_true",  default=False,
                help=""
                )
    parser.add_option("--section-slip-map", "--ssm",
                dest="section_slip_map", type="int",  default=None,
                help="Doc.",
                metavar="WORD_LENGTH"
                )
    parser.add_option("--slip-vs-time", "--svt",
                dest="slip_vs_time", action="store_true",  default=False,
                help="Print event info."
                )
    parser.add_option("--plot-greens-interactions", "--pgi",
                dest="plot_greens_interactions", choices=('both', 'text', 'plot'),  default=False,
                help="",
                metavar="OUTPUT_TYPE"
                )
    parser.add_option("--mags-by-section", "--mbs",
                dest="mags_by_section", action="store_true",  default=False,
                help=""
                )
    parser.add_option("--recurrence-intervals", "--ri",
                dest="recurrence_intervals", choices=('both', 'text', 'plot'),  default=False,
                help="Plots the recurrence intervals for each section in the sim. OUTPUT_TYPE can be text, plot or both.",
                metavar="OUTPUT_TYPE"
                )


    (options, args) = parser.parse_args()

    #the data file
    data_file       = options.data_file
    data_file_type  = options.data_file_type
    fault_trace     = options.fault_trace
    sys_name        = options.sys_name
    debug_mode      = options.debug_mode
    
    
    #these are filters
    if options.section_filter is not None:
        section_filter  = map(int,options.section_filter.split(','))
    else:
        section_filter = options.section_filter
    start_year      = options.start_year
    end_year        = options.end_year
    start_event     = options.start_event
    end_event       = options.end_event
    if options.area_filter is not None:
        if options.area_filter == 'full':
            area_filter = 'full'
        else:
            area_filter_dat  = map(float,options.area_filter.split(','))
            NE_pt = MyVector.LatLonPoint( area_filter_dat[0],area_filter_dat[1])
            SW_pt = MyVector.LatLonPoint( area_filter_dat[2],area_filter_dat[3])
            area_filter = [NE_pt, SW_pt]
    else:
        area_filter = options.area_filter
    movie_only = options.movie_only
    fringes = options.fringes
    if options.look_angles is not None:
        azi, ele = map(float,options.look_angles.split(','))
        look_angles = {'azimuth': azi, 'elevation': ele}
    else:
        look_angles = options.look_angles
    
    
    #these are actions
    plot_faults                         = options.plot_faults
    single_event_animation              = options.single_event_animation
    export_eqsim_geometry               = options.export_eqsim_geometry
    event_animation                     = options.event_animation
    rupture_map                         = options.rupture_map
    frequency_magnitude                 = options.frequency_magnitude
    magnitude_rupture_area              = options.magnitude_rupture_area
    magnitude_average_slip              = options.magnitude_average_slip
    average_slip_rupture_length         = options.average_slip_rupture_length
    export_eqsim_events                 = options.export_eqsim_events    
    eqsim_events_min_slip_map_mag       = float(options.eqsim_events_min_slip_map_mag)
    plot_format                         = options.plot_format
    element_CFF                         = options.element_CFF
    element_CFF_map                     = options.element_CFF_map
    element_CFF_block_entropy           = options.element_CFF_block_entropy
    element_CFF_delta_block_entropy     = options.element_CFF_delta_block_entropy
    element_CFF_E_Hmu                   = options.element_CFF_E_Hmu
    event_displacement_map              = options.event_displacement_map
    event_displacement_map_all          = options.event_displacement_map_all
    event_displacement_map_animation    = options.event_displacement_map_animation
    print_event_info                    = options.print_event_info
    event_rupture_map                   = options.event_rupture_map
    event_time_series                   = options.event_time_series
    event_time_series_marked            = options.event_time_series_marked
    plot_fault_map                      = options.plot_fault_map
    slip_in_slip_out                    = options.slip_in_slip_out
    section_slip_map                    = options.section_slip_map
    slip_vs_time                        = options.slip_vs_time
    plot_greens_interactions            = options.plot_greens_interactions
    mags_by_section                     = options.mags_by_section
    event_fault_map                     = options.event_fault_map
    recurrence_intervals                = options.recurrence_intervals

    do_event_actions = False
    do_IT_actions = False
    
    if options.options_file is not None:
        #parse the options file
        pass
    
    if rupture_map                          is not False or \
            frequency_magnitude             is not False or \
            magnitude_rupture_area          is not False or \
            magnitude_average_slip          is not False or \
            average_slip_rupture_length     is not False or \
            export_eqsim_events             is not False or \
            event_animation                 is not False or \
            print_event_info                is not False or \
            slip_in_slip_out                is not False or \
            recurrence_intervals            is not False:
        do_event_actions = True
    
    
    if element_CFF                          is not False or \
            element_CFF_map                 is not False or \
            element_CFF_block_entropy       is not None  or \
            element_CFF_delta_block_entropy is not None  or \
            element_CFF_E_Hmu               is not None:
        do_IT_actions = True
    
    #print do_event_actions
    
    #if data_file is not None:
    #    try:
    #        os.chdir(os.path.dirname(data_file))
    #    except OSError:
    #        print "error: The data file can not be found"
    #        data_file = None
    
    if fault_trace is not None:
        try:
            os.chdir(os.path.dirname(fault_trace))
        except OSError:
            print "error: The fault trace file can not be found"
            fault_trace = None
    
    if fault_trace is not None:
        if sys_name is None:
            #we need a sys name. try to get it from the name of the data file.
            data_file_split = os.path.basename(fault_trace).split('.')
            sys_name = data_file_split[0]

        export_eqsim_geometry = True
        
    
    if data_file is not None:
        if sys_name is None:
            #we need a sys name. try to get it from the name of the data file.
            data_file_split = os.path.basename(data_file).split('.')
            tmp_sys_name = data_file_split[0]
            
            if data_file_type == 'eqsimg':
                try:
                    end_index = tmp_sys_name.lower().rindex('geometry') - 1
                except ValueError:
                    end_index = len(tmp_sys_name)
            elif data_file_type == 'eqsime':
                try:
                    end_index = tmp_sys_name.lower().rindex('events') - 1
                except ValueError:
                    end_index = len(tmp_sys_name)
            else:
                end_index = len(tmp_sys_name)
            
            sys_name = tmp_sys_name[0:end_index]
    
    if sys_name is None:
        #prompt the user for a system name
        pass
    
    '''this is the correct way to do this. disabled for testing purposes '''
    #vc_sys = VC.VCSys()
    #vc_sys.name = sys_name
    
    #try:
    #    vc_sys.geometry = cPickle.load(open('%s_cache/%s_Geometry.pkl'%sys_name, 'rb'))
    
    #print sys_name
    
    create_sys = True
    #if not debug_mode:
    #    try:
    #        vc_sys = cPickle.load(open('%s_cache/%s.pkl'%(sys_name,sys_name), 'rb'))
    #        create_sys = False
    #    except:
    #        create_sys = True
    
    if create_sys:
        #if not os.path.exists('%s_cache'%sys_name):
        #    os.makedirs('%s_cache'%sys_name)
        vc_sys = VC.VCSys(sys_name, data_file)
        '''
        if data_file is not None:
            if data_file_type == 'h5':
                vc_sys.loadHdf5Geometry(data_file)
                if not debug_mode:
                    vc_sys.eventYearMappingHdf5(data_file)
            if not debug_mode:
                vc_sys.layerSections()
            #disable for testing purposes
            if not debug_mode:
                cPickle.dump(vc_sys, open('%s_cache/%s.pkl'%(sys_name,sys_name), 'wb'))
        '''

    vc_sys.section_filter = section_filter

    if event_displacement_map is not None:
        if area_filter is None:
            vc_sys.eventDisplacementMap(event_displacement_map, plot_area='event', fringes=fringes, look_angles=look_angles)
        else:
            vc_sys.eventDisplacementMap(event_displacement_map, plot_area=area_filter, fringes=fringes, look_angles=look_angles)

    if event_displacement_map_all:
        vc_sys.eventDisplacementMapAll()

    if event_displacement_map_animation:
        vc_sys.eventDisplacementMapAnimation(start_event, end_event, movie_only=movie_only, fringes=fringes, look_angles=look_angles)

    if event_rupture_map is not None:
        vc_sys.eventRuptureMap(event_rupture_map)

    if event_time_series:
        if start_year is not None and start_event == 0:
            start_event = vc_sys.eventForYear(float(start_year))
        
        if end_year is not None and end_event is None:
            end_event = vc_sys.eventForYear(float(end_year))

        
        vc_sys.eventTimeSeries(start_event=start_event, end_event=end_event)

    if event_time_series_marked is not None:
        vc_sys.eventTimeSeries(marked_event=event_time_series_marked)
    
    if plot_fault_map:
        vc_sys.plotFaultMap()

    if event_fault_map is not None:
        vc_sys.plotFaultMap(evid=event_fault_map)

    #if slip_in_slip_out:
    #    if start_year is not None and start_event == 0:
    #        start_event = vc_sys.eventForYear(float(start_year))
    #
    #    if end_year is not None and end_event is None:
    #        end_event = vc_sys.eventForYear(float(end_year))
    #
    #    vc_sys.slipInSlipOut(start_event, end_event)

    if section_slip_map is not None:
        if start_year is not None and start_event == 0:
            start_event = vc_sys.eventForYear(float(start_year))
        
        if end_year is not None and end_event is None:
            end_event = vc_sys.eventForYear(float(end_year))

        vc_sys.sectionSlipMap(section_slip_map, start_event, end_event )

    if slip_vs_time:
        if start_year is not None and start_event == 0:
            start_event = vc_sys.eventForYear(float(start_year))
        
        if end_year is not None and end_event is None:
            end_event = vc_sys.eventForYear(float(end_year))
        
        vc_sys.slipVsTime(start_event, end_event )

    if plot_greens_interactions != False:
        vc_sys.plotGreensInteractions(plot_greens_interactions)

    if mags_by_section:
        if start_year is not None and start_event == 0:
            start_event = vc_sys.eventForYear(float(start_year))
        
        if end_year is not None and end_event is None:
            end_event = vc_sys.eventForYear(float(end_year))
        
        vc_sys.magsBySection(start_event, end_event )
        
    '''the following is just for testing purposes
    vc_sys = VC.VCSys()
    vc_sys.name = sys_name
    vc_sys.loadHdf5Geometry(data_file)
    #vc_sys.proximityTable()
    vc_sys.layerSections()
    #vc_sys.eventYearMappingHdf5(data_file)
    #vc_sys.layerSections()
    #cPickle.dump(vc_sys.geometry, open('%s_Geometry-info.pkl'%sys_name, 'wb'))
    #vc_sys.geometry = cPickle.load(open('%s_Geometry-info.pkl'%sys_name, 'rb'))
    #event_info = cPickle.load(open('%s_Event-info.pkl'%sys_name, 'rb'))
    #vc_sys.events.years = event_info[0]
    #vc_sys.events.event_ids = event_info[1]
    '''
    
    
    
    if do_event_actions:
        if start_year is not None and start_event == 0:
            start_event = vc_sys.eventForYear(float(start_year))
        
        if end_year is not None and end_event is None:
            end_event = vc_sys.eventForYear(float(end_year))
        
        if start_event is None:
            #there was a problem with the start year
            print "error: The given start year is not in the event set."
        
        vc_sys.output_format = plot_format

        #### Rupture Map
        vc_sys.rupture_map                          = rupture_map
        
        #### Event Animation
        vc_sys.event_animation                      = event_animation
        
        #### Frequency Magnitude
        vc_sys.frequency_magnitude                  = frequency_magnitude

        #### Magnitude Rupture Area
        vc_sys.magnitude_rupture_area               = magnitude_rupture_area

        #### Magnitude Average Slip
        vc_sys.magnitude_average_slip               = magnitude_average_slip

        #### Average Slip Rupture Length
        vc_sys.average_slip_rupture_length          = average_slip_rupture_length

        #### Export EqSim Events
        vc_sys.export_eqsim_events_min_slip_map_mag = eqsim_events_min_slip_map_mag
        vc_sys.export_eqsim_events                  = export_eqsim_events
        
        #### Print Event Info
        vc_sys.print_event_info                     = print_event_info

        #### Slip In - Slip Out
        vc_sys.slip_in_slip_out                     = slip_in_slip_out

        #### Recurrence intervals
        vc_sys.recurrence_intervals                 = recurrence_intervals

        #cProfile.runctx('vc_sys.eventActions(data_file, start_event, end_event)', globals(), locals(), 'prof')

        vc_sys.eventActions(start_event, end_event)
    
    if do_IT_actions:
        if start_year is not None and start_event == 0:
            start_event = vc_sys.eventForYear(float(start_year))
        
        if end_year is not None and end_event is None:
            end_event = vc_sys.eventForYear(float(end_year))
        
        if start_event is None:
            #there was a problem with the start year
            print "error: The given start year is not in the event set."
        
        vc_sys.output_format = plot_format
        
        #### Plot Element CFF
        vc_sys.element_CFF                      = element_CFF
        
        #### Plot Element CFF maps
        vc_sys.element_CFF_map                  = element_CFF_map
        
        #### Plot Element CFF block entropies
        vc_sys.element_CFF_block_entropy        = element_CFF_block_entropy
        vc_sys.element_CFF_delta_block_entropy  = element_CFF_delta_block_entropy
        
        #### Plot Element excess entropy and entropy rate
        vc_sys.element_CFF_E_Hmu                = element_CFF_E_Hmu
        
        vc_sys.infoTheoryActions(start_event, end_event)

    if fault_trace is not None:
        vc_sys.geometryFromTraceFile(fault_trace)
    
    if export_eqsim_geometry:
        vc_sys.exportEqsimGeometry('%s_Geometry.dat'%vc_sys.name, '%s_Friction.dat'%vc_sys.name)
    
    if plot_faults:
        vc_sys.exportKMLGeometry('%s_Faults.kml'%vc_sys.name)
    
if __name__ == "__main__": 
    sys.exit(main())
