#!/usr/bin/env python
import tables
import Queue
import multiprocessing
from operator import itemgetter
import time
import itertools

class DataProcessor(multiprocessing.Process):
    def __init__(self, sim_file_path, work_queue, result_queue):
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        
        self.sim_file_path = sim_file_path
        
        super(DataProcessor, self).__init__()
    
    def run(self):
        self.sim_file = tables.open_file('demo_big.h5', 'r')

        sweep_table = self.sim_file.get_node('/event_sweep_table')
        event_table = self.sim_file.get_node('/event_table')

        while not self.kill_received:
 
            # get a task
            try:
                event_range = self.work_queue.get_nowait()
            except Queue.Empty:
                break
            
            #sections = self.vc_sys.geometry.sectionsSortedByID(section_list)
        
            print 'processing events {} - {}'.format(*event_range)
            results = []
            for evnum in range(*event_range):
                areas = {}
                total_slip = 0.0
                slip_records = 0
                for sweep in sweep_table[event_table[evnum]['start_sweep_rec']:event_table[evnum]['end_sweep_rec']]:
                    areas[sweep['block_id']] = sweep['area']
                    total_slip += sweep['slip']
                    slip_records += 1
                #print slip_records
                results.append((evnum, total_slip, slip_records, sum(areas.values())))
            
                #event_elements[i] = set(events.get_event_elements(i))
            
            # the actual processing
            #layered_sections = {}
            #print '    layering %i sections '%(len(sections))
            #for n, s in enumerate(sections):
            #    sid = s.sid
                #print sid
            #    curr_section = []
            #    for eid in s.selement_ids:
            #        if self.vc_sys.geometry.elements[eid].onTrace():
                        #print s.elementsAlongDip(eid)
            #            curr_section.append( s.elementsAlongDip(eid) )
                
            #    curr_section_np = np.array(curr_section)
                
            #    layered_sections[s.sid] = curr_section_np.T

            # store the result
            self.result_queue.put(results)
        
        self.sim_file.close()


class VCSimData(object):
    def __init__(self, file_path=None):
        self.file = None
        self.file_path = file_path
        self.do_event_area = False
        self.do_event_average_slip = False
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

        if self.do_event_area or self.do_event_average_slip:
            self.calculate_additional_data()

    def calculate_additional_data(self):
    # get the new data
        print 'getting new data'
        start_time = time.time()
        num_processes = multiprocessing.cpu_count()
        total_events = self.file.root.event_table.nrows
        #close the current file
        self.file.close()
        
        #break the work up
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
            worker = DataProcessor(self.file_path, work_queue, result_queue)
            worker.start()

        # collect the results off the queue
        results = []
        for i in range(num_processes):
            results.append(result_queue.get())
        data_process_results_sorted = sorted(list( itertools.chain(*results) ), key=itemgetter(0))
        
        print 'Done! {} seconds'.format(time.time() - start_time)
        #print data_process_results_sorted

    #create the new table
        self.file = tables.open_file(self.file_path, 'a')
        table = self.file.root.event_table
        
        # Get a description of table in dictionary format
        descr = table.description._v_colObjects
        descr2 = descr.copy()
        
        # Add a column to description
        if self.do_event_area:
            descr2['event_area'] = tables.Float64Col(dflt=0.0)
        
        if self.do_event_average_slip:
            descr2['event_average_slip'] = tables.Float64Col(dflt=0.0)
        
        # Create a new table with the new description
        table2 = self.file.create_table('/', 'tmp', descr2, 'Event Table')
        
        # Copy the user attributes
        table.attrs._f_copy(table2)
        
        # Fill the rows of new table with default values
        for i in xrange(table.nrows):
            table2.row.append()
        
        # Flush the rows to disk
        table2.flush()
        
        # Copy the columns of source table to destination
        for col in descr:
            getattr(table2.cols, col)[:] = getattr(table.cols, col)[:]

        # Fill the new column
        if self.do_event_area:
            table2.cols.event_area[:] = [ x[3] for x in data_process_results_sorted ]
        
        if self.do_event_average_slip:
            table2.cols.event_average_slip[:] = [ x[1]/float(x[2]) for x in data_process_results_sorted ]
        	
        # Remove the original table
        table.remove()
        
        # Move table2 to table
        table2.move('/','event_table')
        	
        # Print the new table
        #print "Contents of the table with column added:", self.file.root.event_table[:]
        
        # Finally, close the file
        self.file.close()

    #open the file with the new table
        self.file = tables.open_file(self.file_path)

        #self.file.root.event_table.col('event_num')
'''
# All access to the file goes through a single instance of this class.
# It contains several queues that are used to communicate with other
# processes.
# The read_queue is used for requests to read data from the HDF5 file.
# A list of result_queues is used to send data back to client processes.
# The write_queue is used for requests to modify the HDF5 file.
# One end of a pipe (shutdown) is used to signal the process to terminate.
class FileAccess(multiprocessing.Process):

    def __init__(self, h5_path, node_name, read_queue, result_queues, write_queue, shutdown):
        self.h5_path = h5_path
        self.node_name = node_name
        self.read_queue = read_queue
        self.result_queues = result_queues
        self.write_queue = write_queue
        self.shutdown = shutdown
        self.block_period = .01
        super(FileAccess, self).__init__()

    def run(self):
        self.h5_file = tables.open_file(self.h5_path, 'r+')
        self.table = self.h5_file.get_node('/{}'.format(self.node_name))
        another_loop = True
        while another_loop:

            # Check if the process has received the shutdown signal.
            if self.shutdown.poll():
                another_loop = False

            # Check for any data requests in the read_queue.
            try:
                row_num, proc_num = self.read_queue.get(True, self.block_period)
                # look up the appropriate result_queue for this data processor
                # instance
                result_queue = self.result_queues[proc_num]
                print 'processor {0} reading from row {1}'.format(proc_num, row_num)
                result_queue.put(self.read_data(row_num))
                another_loop = True
            except Queue.Empty:
                pass

            # Check for any write requests in the write_queue.
            try:
                row_num, data = self.write_queue.get(True, self.block_period)
                print 'writing row', row_num
                self.write_data(row_num, data)
                another_loop = True
            except Queue.Empty:
                pass

        # close the HDF5 file before shutting down
        self.h5_file.close()

    def read_data(self, row_num):
        return self.table[row_num, :]

    def write_data(self, row_num, data):
        self.table[row_num, :] = data

# This class represents a process that does work by reading and writing to the
# HDF5 file.  It does this by sending requests to the FileAccess class instance
# through its read and write queues.  The data results are sent back through
# the result_queue.
# Its actions are logged to a text file.
class DataProcessor(multiprocessing.Process):

    def __init__(self, read_queue, result_queue, write_queue, proc_num):
        self.read_queue = read_queue
        self.result_queue = result_queue
        self.write_queue = write_queue
        self.proc_num = proc_num
        #self.output_file = output_file
        super(DataProcessor, self).__init__()

    def run(self):
        #self.output_file = open(self.output_file, 'w')
        # read a random row from the file
        row_num = random.randint(0, self.array_size - 1)
        self.read_queue.put((row_num, self.proc_num))
        self.output_file.write(str(row_num) + '\n')
        self.output_file.write(str(self.result_queue.get()) + '\n')

        # modify a random row to equal 11 * (self.proc_num + 1)
        #row_num = random.randint(0, self.array_size - 1)
        #new_data = (numpy.zeros((1, self.array_size), 'i8') +
        #            11 * (self.proc_num + 1))
        #self.write_queue.put((row_num, new_data))

        # pause, then read the modified row
        #time.sleep(0.015)
        #self.read_queue.put((row_num, self.proc_num))
        #self.output_file.write(str(row_num) + '\n')
        #self.output_file.write(str(self.result_queue.get()) + '\n')
        #self.output_file.close()
'''
