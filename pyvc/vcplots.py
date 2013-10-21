#!/usr/bin/env python
from pyvc import *
from pyvc import vcutils
from pyvc import vcexceptions
import matplotlib.pyplot as mplt
import matplotlib.font_manager as mfont
import numpy as np
import math
import multiprocessing
import cPickle
import networkx as nx

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