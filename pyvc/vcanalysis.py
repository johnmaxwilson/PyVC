from pyvc import *
from pyvc import vcutils
from operator import itemgetter
import networkx as nx
from subprocess import call
import cPickle
import sys
import numpy as np
import matplotlib.pyplot as mplt
import itertools
from collections import deque

def cum_prob(sim_file, output_file, event_range=None, section_filter=None, magnitude_filter=None,plot_type=1):
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)

        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)

        event_data = events.get_event_data(['event_number', 'event_year', 'event_magnitude', 'event_range_duration'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)

    intervals = [
                x - event_data['event_year'][n-1]
                for n,x in enumerate(event_data['event_year'])
                if n != 0]
    
    # t vs P(t)
    if plot_type ==1:
        mplt.plot([x for x in sorted(intervals)], [float(n)/float(len(intervals)) for n,x in enumerate(sorted(intervals))])
    # t0 vs P(t0 + dt, t0)
    elif plot_type == 2:
        dt = 100
    
        t0s = []
        pts = []
    
        for t0 in [0.0] + [x for x in sorted(intervals)]:
            intervals = [x - event_data['event_year'][n-1] for n,x in enumerate(event_data['event_year']) if n != 0 and x - event_data['event_year'][n-1] > t0+dt]
            intervals_t0 = [x - event_data['event_year'][n-1] for n,x in enumerate(event_data['event_year']) if n != 0 and x - event_data['event_year'][n-1] > t0]
        
            if len(intervals_t0) != 0:
                t0s.append(t0)
                pts.append(1.0 - float(len(intervals))/float(len(intervals_t0)))
    
        mplt.plot(t0s,pts)
    # t=t0+dt vs P(t,t0)
    elif plot_type == 3:
        for t0 in range(0,175,25):
            ts = []
            P_t_t0 = []
            intervals_t0 = [x - event_data['event_year'][n-1] for n,x in enumerate(event_data['event_year']) if n != 0 and x - event_data['event_year'][n-1] > t0]
            for dt in range(250):
                intervals = [x - event_data['event_year'][n-1] for n,x in enumerate(event_data['event_year']) if n != 0 and x - event_data['event_year'][n-1] > t0+dt]

                if len(intervals_t0) != 0:
                    ts.append(t0+dt)
                    P_t_t0.append(1.0 - float(len(intervals))/float(len(intervals_t0)))

            mplt.plot(ts,P_t_t0)
            
    else:
        sys.exit("Error, choose plot_type = 1, 2, or 3!")
        
    mplt.savefig(output_file,dpi=100)
    sys.stdout.write("\Plot saved to: "+output_file+'\n')
    
    return event_data['event_year']

#-------------------------------------------------------------------------------
def weibull(x_array,beta,tau):
    if len(x_array) < 2:
        sys.exit("Input must be an array")
    else:
        return np.array([1-np.exp( -(x/float(tau))**beta) for x in x_array])
        
#-------------------------------------------------------------------------------
def cond_weibull(x_array,t0,beta,tau,single=False):
    if single==False:
        if len(x_array) < 2:
            sys.exit("Input must be an array")
        else:
            return np.array([1-np.exp( (t0/float(tau))**beta - (x/float(tau))**beta) for x in x_array])
    else:
        return 1-np.exp( (t0/float(tau))**beta - (x_array/float(tau))**beta)

#-------------------------------------------------------------------------------
def cond_weibull_fixed_dt(x_array,dt,beta,tau):
    if len(x_array) < 2:
        sys.exit("Input must be an array")
    else:
        return np.array([1-np.exp( (x/float(tau))**beta - ((x+dt)/float(tau))**beta) for x in x_array])

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

#-------------------------------------------------------------------------------
def eq_info(sim_file, evnums):

     with VCSimData() as sim_data:
     
        # open the simulation data file
        sim_data.open_file(sim_file)

        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)

        requested_data = ['event_number', 'event_year', 'event_magnitude','event_average_slip','event_surface_rupture_length']
        
        event_data = events.get_event_data_from_evnums(evnums, requested_data)
        
        print '\n{0:<10}{1:<10}{2:<10}{3:<10}{4:<10}'.format('num','year','magnitude','slip [m]','rupt.len [km]')

        for i in range(len(evnums)):
            print '{ev_num:<10}{ev_year:<10.2f}{ev_mag:<10.2f}{ev_av_slip:<10.2f}{ev_rup_len:<10.2f}'.format(ev_num=event_data['event_number'][i], ev_year=event_data['event_year'][i], ev_mag=event_data['event_magnitude'][i],ev_av_slip=event_data['event_average_slip'][i],ev_rup_len=event_data['event_surface_rupture_length'][i]/1000.0)
        print "\n"
  



#-------------------------------------------------------------------------------
def detailed_sim_info(sim_file, sortby='event_magnitude', show=15, event_range=None, section_filter=None, magnitude_filter=None,return_evnums=False,min_mag=0.0):

     with VCSimData() as sim_data:
        evnums = []     
     
        # open the simulation data file
        sim_data.open_file(sim_file)

        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)

        event_data = events.get_event_data(['event_number', 'event_year', 'event_magnitude', 'event_range_duration','event_average_slip','event_surface_rupture_length'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        
        print '{0:<10}{1:<10}{2:<10}{3:<10}{4:<10}'.format('num','year','magnitude','slip [m]','rupt.len [km]')
        if sortby == 'event_elements':
            sorted_data = [i[0] for i in sorted(enumerate(event_data[sortby]), lambda a,b: cmp(len(b[1]),len(a[1])), reverse=True)][0:show]
        else:
            sorted_data = [i[0] for i in sorted(enumerate(event_data[sortby]), key=itemgetter(1), reverse=True)][0:show]
     
        for i in sorted_data:
            if event_data['event_magnitude'][i] > min_mag:
                print '{ev_num:<10}{ev_year:<10.2f}{ev_mag:<10.2f}{ev_av_slip:<10.2f}{ev_rup_len:<10.2f}'.format(ev_num=event_data['event_number'][i], ev_year=event_data['event_year'][i], ev_mag=event_data['event_magnitude'][i],ev_av_slip=event_data['event_average_slip'][i],ev_rup_len=event_data['event_surface_rupture_length'][i]/1000.0)
                if return_evnums:
                    evnums.append(event_data['event_number'][i])

     if return_evnums:
        return evnums



def event_sections(sim_file,evnum):
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)

        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)

        elements = events.get_event_elements(evnum)
        sections = geometry.sections_with_elements(elements)
        
        sys.stdout.write('\nevent {}\n'.format(evnum))    
        for secid in sections:
            sys.stdout.write('{sec:<3}{name:<10}\t'.format(sec=secid,name=geometry.get_section_name(secid)))
    




def graph_events(sim_file, output_file, triggers_only=False, event_range=None, section_filter=None, magnitude_filter=None):
    
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
        if triggers_only:
            event_data = events.get_event_data(['event_trigger', 'event_year', 'event_magnitude', 'event_number'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        else:
            event_data = events.get_event_data(['event_elements', 'event_year', 'event_magnitude', 'event_number'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        
        # initilize a graph
        G = nx.DiGraph(name='Event graph for {}'.format(sim_file), sim_file=sim_file, event_range=None, section_filter=None, magnitude_filter=None)
        
        sys.stdout.write('{} events : {} years\n'.format(len(event_data['event_year']),event_data['event_year'][-1] - event_data['event_year'][0] ))
        sys.stdout.flush()
        
        # add edges and nodes to the graph for each event
        if triggers_only:
            ev_elements = [[x] for x in event_data['event_trigger']]
        else:
            ev_elements = event_data['event_elements']
        for i, ev_eles in enumerate(ev_elements):
            
            if i%round(float(len(event_data['event_year']))/100.0) == 0:
                sys.stdout.write('\r event {} of {}'.format(i, len(event_data['event_year'])))
                sys.stdout.flush()
            for this_sid in geometry.sections_with_elements(ev_eles):
                if i < len(ev_elements) - 1:
                    for next_sid in geometry.sections_with_elements(ev_elements[i+1]):
                        duration = event_data['event_year'][i+1] - event_data['event_year'][i]
                        try:
                            G[this_sid][next_sid]['weight'] += 1
                            G[this_sid][next_sid]['duration'].append(duration)
                        except KeyError:
                            G.add_edge(this_sid, next_sid, weight=1, duration=[duration])
                        G.node[next_sid]['type'] = 'section'
                    G.node[this_sid]['magnitude'] = event_data['event_magnitude'][i]
                    G.node[this_sid]['number'] = event_data['event_number'][i]
                    G.node[this_sid]['type'] = 'section'
    
        # add the duration mean and standard deviation
        for i in G:
            for j in G[i]:
                G[i][j]['duration_mean'] = np.mean(G[i][j]['duration'])
                G[i][j]['duration_std'] = np.std(G[i][j]['duration'])

        # save the graph
        sys.stdout.write('\nSaving graph ')
        sys.stdout.flush()
        cPickle.dump(G, open(output_file, 'wb'))

def analyze_event_sequence_graph(graph_file):
    G = cPickle.load(open(graph_file, 'rb'))
    
    sequences_by_degree = {}
    for n in nx.nodes_iter(G):
        if G.node[n]['type'] == 'section':
            sequences_by_degree[n] = G.degree(n)

    sorted_seq = sorted(sequences_by_degree.iteritems(), key=itemgetter(0))

    print sorted_seq

    # plot parameters
    imw = 1024.0 # the full image width
    imh = 1024.0
    lm = 40.0
    rm = 50.0
    tm = 50.0
    bm = 50.0
    res = 72.0
    
    imwi = imw/res
    imhi = imh/res
    fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
    ph = imh - tm - bm # the height for both matricies
    pw = imw - lm - rm
    ax = fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))

    ax.plot(range(len(sorted_seq)),[x[1] for x in sorted_seq])

    print [x for x in G.edges(sorted_seq[0][0], data=True)]
    print sorted_seq[0][0]


def graph_event_sequences(sim_file, output_file, sequence_length=5, event_range=None, section_filter=None, magnitude_filter=None):
    
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
        event_data = events.get_event_data(['event_elements', 'event_year', 'event_magnitude', 'event_number', 'event_trigger'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        
        # initilize a graph
        G = nx.DiGraph(name='Event sequence graph for {}'.format(sim_file), sim_file=sim_file, event_range=None, section_filter=None, magnitude_filter=None)
        
        sys.stdout.write('{} events : {} years\n'.format(len(event_data['event_year']),event_data['event_year'][-1] - event_data['event_year'][0] ))
        sys.stdout.flush()
        
        # add edges and nodes to the graph for each event
        current_sequence = deque()
        for i, event_trigger in enumerate(event_data['event_trigger']):
            trigger_sid = geometry.sections_with_elements([event_trigger])[0]
            if i > sequence_length - 1:
                this_sequence_label = '->'.join([str(x) for x in current_sequence])
                #print this_sequence_label
                if i < len(event_data['event_trigger']) - 1:
                    for next_sid in geometry.sections_with_elements(event_data['event_elements'][i+1]):
                        duration = event_data['event_year'][i+1] - event_data['event_year'][i]
                        magnitude = event_data['event_magnitude'][i+1]
                        try:
                            G[this_sequence_label][next_sid]['weight'] += 1
                            G[this_sequence_label][next_sid]['duration'].append(duration)
                            G[this_sequence_label][next_sid]['magnitude'].append(magnitude)
                        except KeyError:
                            G.add_edge(this_sequence_label, next_sid, weight=1, duration=[duration], magnitude=[magnitude])
                
                        G.node[next_sid]['type'] = 'section'
                        G.node[next_sid]['bipartite'] = 0
            
                    G.node[this_sequence_label]['type'] = 'sequence'
                    G.node[this_sequence_label]['bipartite'] = 1
        
                current_sequence.popleft()
            current_sequence.append(trigger_sid)
            
            '''
            if i%sequence_length == 0:
                if len(current_sequence) == 0:
                    current_sequence.append(this_sid)
                else:
                    this_sequence_label = '-'.join([str(x) for x in current_sequence])
                    if last_sequence_label is not None:
                        try:
                            G[last_sequence_label][this_sid]['weight'] += 1
                        except KeyError:
                            G.add_edge(last_sequence_label, this_sid, weight=1)
                    last_sequence_label = this_sequence_label
                    current_sequence = [this_sid]
            else:
                current_sequence.append(this_sid)
            '''
        
        
            '''
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
                '''
        
        '''
        # add the duration mean and standard deviation
        for i in G:
            for j in G[i]:
                G[i][j]['duration_mean'] = np.mean(G[i][j]['duration'])
                G[i][j]['duration_std'] = np.std(G[i][j]['duration'])
        '''
        # save the graph
        sys.stdout.write('\nSaving graph ')
        sys.stdout.flush()
        cPickle.dump(G, open(output_file, 'wb'))

def generate_event_sequence(graph_file, start_sid, length=100, runs=1):
    G = cPickle.load(open(graph_file, 'rb'))
    
    matrix, pos_sid = nx.attr_matrix(G, edge_attr='weight', normalized=True)
    sid_pos = {sid: position for (position, sid) in enumerate(pos_sid)}
    
    duration_mean_matrix, pos_sid_mean = nx.attr_matrix(G, edge_attr='duration_mean')
    
    raw_output = np.empty((length,runs))
    output = np.empty(length)
    time = np.empty(length)
    
    current_run = 0
    while current_run < runs:
        current_step = 0
        current_time = 0.0
        current_node = start_sid
        while current_step < length:
            
            raw_output[current_step, current_run] = current_node
            time[current_step] = current_time
            
            out_probs = np.cumsum(matrix[sid_pos[current_node]], axis=1)
            
            #for out_node in G[current_node]:
            #    print out_node, matrix[sid_pos[current_node], sid_pos[out_node]], out_probs[0,sid_pos[out_node]]
            
            choice = np.random.random_sample()
            
            choice_index = np.argwhere(out_probs<choice)
            
            try:
                next_node = pos_sid[choice_index[-1,0,-1]+1]
            except IndexError:
                next_node = pos_sid[0]
            
            current_time += duration_mean_matrix[sid_pos[current_node], sid_pos[next_node]]
            
            current_node = next_node
            
            #try:
            #    print choice, choice_index[-1,0,-1], pos_sid[choice_index[-1,0,-1]+1]
            #except IndexError:
            #    print 0, pos_sid[0]
            
            #print choice, choice_index[-1,0,-1], pos_sid[choice_index[-1,0,-1]+1]
            
            current_step += 1
        current_run += 1


    print raw_output.shape

    for index in range(length):
        output[index] = np.mean(raw_output[index])

    xs = []
    ys = []

    for index in range(length):
        if index < length - 1:
            xs.append(output[index])
            ys.append(output[index+1])


    # plot parameters
    imw = 1024.0 # the full image width
    imh = 1024.0
    lm = 40.0
    rm = 50.0
    tm = 50.0
    bm = 50.0
    res = 72.0
    
    imwi = imw/res
    imhi = imh/res
    fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
    ph = imh - tm - bm # the height for both matricies
    pw = imw - lm - rm
    ax = fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))

    #ax.plot(time, output)
    #ax.set_ylim((1, max(pos_sid)))

    ax.scatter(xs,ys)
    ax.set_ylim((0.5, max(pos_sid)+0.5))
    ax.set_xlim((0.5, max(pos_sid)+0.5))


def find_event_sequence_r(sid, matrix, pos_sid, sid_pos, depth, results, stack, top):
    indices =  (np.argsort(matrix[sid_pos[sid], :]).T)[::-1][0:top]
    
    depth -= 1
    
    stack.append(sid)
    
    if depth >= 0:
        for i in indices:
            
            find_event_sequence_r( pos_sid[i[0,0]], matrix, pos_sid, sid_pos, depth, results, stack, top)
        stack.pop()
    else:
        for i in stack:
            results.append(i)
        stack.pop()

def sequence_probability(sequence, matrix, sid_pos):
    
    ret = 1
    
    for i in range(sequence.size):
        try:
            ret *= matrix[sid_pos[sequence[i]], sid_pos[sequence[i+1]]]
        except IndexError:
            pass

    return ret

def find_event_sequence(graph_file, start_sid, length, top=3):
    G = cPickle.load(open(graph_file, 'rb'))
    
    matrix, pos_sid = nx.attr_matrix(G, edge_attr='weight', normalized=True)
    
    sid_pos = {sid: position for (position, sid) in enumerate(pos_sid)}
    
    results = []
    find_event_sequence_r(start_sid, matrix, pos_sid, sid_pos, length, results, [], top)
    
    _results = np.reshape(np.array(results), (-1, length+1))
    
    ret_unsorted = [{'sequence':_results[i], 'probability':sequence_probability(_results[i], matrix, sid_pos)}
                    for i in range(_results.shape[0])]
    
    ret_sorted = [x for x in sorted(ret_unsorted, key=lambda x: x['probability'], reverse=True)]
    
    return ret_sorted
    #for i in range(_results.shape[0]):
    #    print _results[i], sequence_probability(_results[i], matrix, sid_pos)
    
    #print _results[0]
    #print len(results), _results.shape, _results.size
    #indices =  (np.argsort(matrix[start_sid, :]).T)[::-1][0:3]
    #print indices[::-1]
    #for i in indices:
    #    print i[0,0]
    
    #for i in itertools.permutations(order,length):
    #    print i
    #print node_map
    '''
    my_matrix = np.zeros((len(G), len(G)))
    
    node_map = {node: key for (key, node) in enumerate(G)}
    #for i, node in enumerate(G):
    
    for i, node in enumerate(G):
        for sid, info in G[node].iteritems():
            j = node_map[sid]
            my_matrix[i,j] = info['weight']
    
    n1 = 10
    n2 = 10
    
    total_weights = 0
    for sid, info in G[n1].iteritems():
        total_weights += info['weight']
    
    print my_matrix[n1, n2], matrix[n1, n2], total_weights
    '''
    '''
    for node in G:
        print node,
        for neighbor in G[node]:
            print neighbor,
        print
    '''
    
    
    #print '11,10', matrix[11, 10]
    #print '10,11', matrix[10, 11]
    #print order
    
    #it = np.nditer(matrix[start_sid, 0:20], flags=['c_index'])
    #while not it.finished:
    #    print it.index, order[it.index], it[0]
    #    it.iternext()
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
    
    shear_ax.imshow(matrix.T, interpolation='none')
    #shear_ax.axis('tight')
    #shear_ax.set_ylim((15.5, 0.5))
    #shear_ax.set_xlim((0.5, 15.5))
    '''
    
    
    '''
    fig.savefig('local/graph_matrix.png', format='png')
    '''
    '''
    total_weights = 0
    for sid, info in G[start_sid].iteritems():
        total_weights += info['weight']
    
    for sid, info in G[start_sid].iteritems():
        print sid, float(info['weight'])/float(total_weights), np.mean(info['duration']), np.std(info['duration']), start_sid
    '''

