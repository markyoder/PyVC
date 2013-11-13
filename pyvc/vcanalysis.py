from pyvc import *
from pyvc import vcutils
from operator import itemgetter
import networkx as nx
from subprocess import call
import cPickle
import sys
import numpy as np
import matplotlib.pyplot as mplt

#-------------------------------------------------------------------------------
# Prints out various information about a simulation.
#-------------------------------------------------------------------------------
def sim_info(sim_file, sortby='event_magnitude', show=50, event_range=None, section_filter=None, magnitude_filter=None):
     with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)

        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)

        event_data = events.get_event_data(['event_number', 'event_year', 'event_magnitude', 'event_range_duration'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        
        print '{0:<10}{1:<10}{2:<10}'.format('num','year','magnitude')
        if sortby == 'event_elements':
            sorted_data = [i[0] for i in sorted(enumerate(event_data[sortby]), lambda a,b: cmp(len(b[1]),len(a[1])), reverse=True)][0:show]
        else:
            sorted_data = [i[0] for i in sorted(enumerate(event_data[sortby]), key=itemgetter(1), reverse=True)][0:show]
     
        for i in sorted_data:
            print '{ev_num:<10}{ev_year:<10.2f}{ev_mag:<10.2f}'.format(ev_num=event_data['event_number'][i], ev_year=event_data['event_year'][i], ev_mag=event_data['event_magnitude'][i])

def graph_events(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None):
    
    sys.stdout.write('Initializing graph :: ')
    sys.stdout.flush()
    
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)

        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)
        
        # get the data
        event_data = events.get_event_data(['event_elements', 'event_year', 'event_magnitude', 'event_number'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        
        # initilize a graph
        G = nx.DiGraph(sim_file=sim_file, event_range=None, section_filter=None, magnitude_filter=None)
        
        sys.stdout.write('{} events : {} years\n'.format(len(event_data['event_year']),event_data['event_year'][-1] - event_data['event_year'][0] ))
        sys.stdout.flush()
        
        # add edges and nodes to the graph for each event
        for i, ev_eles in enumerate(event_data['event_elements']):
            if i%round(float(len(event_data['event_year']))/100.0) == 0:
                sys.stdout.write('\r event {} of {}'.format(i, len(event_data['event_year'])))
                sys.stdout.flush()
            for this_sid in geometry.sections_with_elements(ev_eles):
                try:
                    for next_sid in geometry.sections_with_elements(event_data['event_elements'][i+1]):
                        duration = event_data['event_year'][i+1] - event_data['event_year'][i]
                        try:
                            G[this_sid][next_sid]['weight'] += 1
                            G[this_sid][next_sid]['duration'].append(duration)
                        except KeyError:
                            G.add_edge(this_sid, next_sid, weight=1, duration=[duration])
                        G.node[this_sid]['magnitude'] = event_data['event_magnitude'][i]
                        G.node[this_sid]['number'] = event_data['event_number'][i]
                except IndexError:
                    pass
    
        # save the graph
        sys.stdout.write('\nSaving graph ')
        sys.stdout.flush()
        cPickle.dump(G, open(output_file, 'wb'))

def event_sequence(graph_file, start_sid, length):
    G = cPickle.load(open(graph_file, 'rb'))
    
    matrix, order = nx.attr_matrix(G, edge_attr='weight', normalized=False)
    
    print G[10][10]['weight'], matrix[10, 10]
    
    '''
    for node in G:
        print node,
        for neighbor in G[node]:
            print neighbor,
        print
    '''
    '''
    
    # plot parameters
    imw = 1024.0 # the full image width
    imh = 1024.0
    lm = 40.0
    rm = 50.0
    tm = 50.0
    bm = 50.0
    res = 72.0
    cbh = 20.0
    cbs = 40.0
    
    #arial14 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=14)
    #arial12 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=12)
    #arial10 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=10)
    #arial7_light = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=7, weight='light')
    
    imwi = imw/res
    imhi = imh/res
    fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
    ph = imh - tm - bm - cbh - cbs # the height for both matricies
    pw = imw - lm - rm
    shear_ax = fig.add_axes((lm/imw, (bm+cbh+cbs)/imh, pw/imw, ph/imh))
    
    shear_ax.imshow(matrix)
    shear_ax.invert_yaxis()
    shear_ax.axis('tight')
    
    fig.savefig('local/graph_matrix.png', format='png')
    '''
    '''
    total_weights = 0
    for sid, info in G[start_sid].iteritems():
        total_weights += info['weight']
    
    for sid, info in G[start_sid].iteritems():
        print sid, float(info['weight'])/float(total_weights), np.mean(info['duration']), np.std(info['duration']), start_sid
    '''

