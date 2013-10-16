from pyvc import *
from pyvc import vcutils
from operator import itemgetter
import networkx as nx
from subprocess import call
import cPickle

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
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)

        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)
        
        # get the data
        event_data = events.get_event_data(['event_elements'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        
        # initilize a graph
        G = nx.DiGraph(event_range=None, section_filter=None, magnitude_filter=None)
        
        # add edges and nodes to the graph for each event
        for i, ev_eles in enumerate(event_data['event_elements']):
            for this_sid in geometry.sections_with_elements(ev_eles):
                try:
                    for next_sid in geometry.sections_with_elements(event_data['event_elements'][i+1]):
                        try:
                            G[this_sid][next_sid]['weight'] += 1
                        except KeyError:
                            G.add_edge(this_sid, next_sid, weight=1)
                except IndexError:
                    pass
    
        # save the graph
        cPickle.dump(G, open(output_file, 'wb'))
