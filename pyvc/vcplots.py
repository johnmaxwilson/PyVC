#!/usr/bin/env python
import matplotlib.pyplot as mplt
mplt.switch_backend('agg')

from pyvc import *
from pyvc import vcutils
from pyvc import vcplotutils
from pyvc import vcexceptions
from pyvc import vcanalysis

import matplotlib.font_manager as mfont
import matplotlib.colors as mcolor
import matplotlib.colorbar as mcolorbar
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap

import numpy as np
import math
import cPickle
import time
from operator import itemgetter
import os
import sys
import gc
import itertools
import subprocess
import networkx as nx
from networkx.algorithms import bipartite

import quakelib



#-------------------------------------------------------------------------------
# New plot to check median recurrence interval vs mode  
#-------------------------------------------------------------------------------
def recurrence_stats(sim_file, event_range=None, section_filter=None,
    magnitude_filter=None, Nbins=20, title_str=""):

    #---------------------------------------------------------------------------
    # Get the data.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)

        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        
        event_data = events.get_event_data(['event_number', 'event_year', 'event_magnitude', 'event_range_duration'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)


    intervals = np.array([   x - event_data['event_year'][n-1]
                    for n,x in enumerate(event_data['event_year'])
                    if n != 0
                ])
    
    # Plot histogram of recurrence intervals
    recurr_ax = mplt.subplot(111)
    recurr_ax.hist(intervals, Nbins,color='g')
    mplt.xlabel("Recurrence interval [years]")
    mplt.title("Recurrence intervals for "+str(len(intervals))+" M "+magnitude_filter+" "+title_str+" EQs ")
    
    # Generate CDF, find median
    cumulative = {}
    cumulative['x'] = np.sort(intervals)
    cumulative['y'] = np.arange(float(intervals.size))/float(intervals.size)
    
    # Find median
    median_index = (np.abs(np.array(cumulative['y'])-0.5)).argmin()
    median_recurr = cumulative['x'][median_index]
    mean_recurrence = np.mean(intervals)
    
    sys.stdout.write("\nMedian recurrence: "+str(median_recurr)+" years")
    sys.stdout.write("\nMean recurrence: "+str(mean_recurrence)+" years\n")

    # Draw verticle lines denoting the mean and median
    recurr_ax.axvline(x=median_recurr,ymax=mplt.ylim()[1], 
            label="median = {:5.2f}".format(median_recurr), color='r')
    recurr_ax.axvline(x=median_recurr,ymax=mplt.ylim()[1],
            label="mean = {:5.2f}".format(mean_recurrence), color='k')
            
    mplt.legend(loc='best')
    
    mplt.savefig('local/test.pdf', dpi=70)




#-------------------------------------------------------------------------------
# Save all forecast plots in one column of many subplots  
#-------------------------------------------------------------------------------
def plot_forecast(sim_file, output_file=None, event_graph_file=None, 
    event_sequence_graph_file=None, event_range=None, section_filter=None,
    magnitude_filter=None, padding=0.08, fixed_dt=30.0):

    #---------------------------------------------------------------------------
    # Plot parameters.
    #---------------------------------------------------------------------------
    imw = 816.0 # The full image width.
    # Image height is set based on the sub plots.
    sph = 200.0 # The sub-plot height
    # Sub-plot width is based on the specific sub-plot.
    lm = 50.0
    rm = 15.0
    tm = 40.0
    bm = 20.0
    sphs = 40.0 # The sub-plot horizontal spacing.
    spvs = 40.0 # The sub-plot vertical spacing.
    
    mmw = 500.0
    
    res = 72.0
    map_res = 'i'
    map_proj = 'cyl'
    
    titlefont1 = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=12)
    titlefont2 = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=14, weight='bold')
    sectionkeyfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=7)
    ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
    framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
    legendfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
    smtitlefont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9, weight='bold')
    cbticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
    
    matrix_cmap = mplt.get_cmap('hot_r')
    sequence_cmap = mplt.get_cmap('autumn')
    
    water_color             = '#bed5ff'
    land_color              = '#ffffff'
    seq_land_color          = '#ffffff'
    boundary_color          = '#000000',
    coastline_color         = '#9a9a9a'
    country_color           = '#9a9a9a'
    state_color             = '#9a9a9a'
    fault_color             = '#000000'
    alt_fault_color         = '#737373'
    selected_fault_color    = '#FFFFFF'
    map_tick_color          = '#000000'
    map_frame_color         = '#000000'
    grid_color              = '#000000'
    cb_fontcolor            = '#000000'

    boundary_width          = 1.0
    coastline_width         = 1.0
    country_width           = 1.0
    state_width             = 1.0
    fault_width             = 0.5
    forecast_fault_width    = 6.0
    seq_fault_width_max     = 6.0
    seq_fault_width_min     = 3.0
    map_frame_width         = 1.0
    grid_width              = 0.5
    num_grid_lines          = 5
    
    sp_line_color           = '#000000'
    sp_line_colormap        = sequence_cmap
    sp_line_width           = 2.0
    
    t0_dt_main_line_color   = '#000000'
    t0_dt_sub_line_color    = '#737373'
    t0_dt_main_line_width   = 2.0
    t0_dt_sub_line_width    = 1.0
    t0_dt_range_color       = sequence_cmap(0.99)
    
    prob_cbh                = 20.0
    prob_cbs                = 40.0
    
    forcast_title_line_spacing = 20.0
    
    legend_loc='best'

    #---------------------------------------------------------------------------
    # Get the data.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)

        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        
        event_data = events.get_event_data(['event_number', 'event_year', 'event_magnitude', 'event_range_duration'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
    
        geometry = VCGeometry(sim_data)
        
        min_lat = geometry.min_lat
        max_lat = geometry.max_lat
        min_lon = geometry.min_lon
        max_lon = geometry.max_lon
        base_lat = geometry.base_lat
        base_lon = geometry.base_lon

        fault_traces = geometry.get_fault_traces()

        section_names = {sid:geometry.get_section_name(sid) for sid in fault_traces.iterkeys()}

    intervals = np.array([   x - event_data['event_year'][n-1]
                    for n,x in enumerate(event_data['event_year'])
                    if n != 0
                ])
    
            

    
    # Calculate the lat-lon range based on the min-max and the padding
    lon_range = max_lon - min_lon
    lat_range = max_lat - min_lat
    max_range = max((lon_range, lat_range))
    min_lon = min_lon - lon_range*padding
    min_lat = min_lat - lat_range*padding
    max_lon = max_lon + lon_range*padding
    max_lat = max_lat + lat_range*padding
    
    # A conversion instance for doing the lat-lon to x-y conversions
    convert = quakelib.Conversion(base_lat, base_lon)
    
    # Convert the fault traces to lat-lon
    fault_traces_latlon = {}
    for secid in fault_traces.iterkeys():
         fault_traces_latlon[secid] = zip(*[(lambda y: (y.lat(),y.lon()))(convert.convert2LatLon(quakelib.Vec3(x[0], x[1], x[2]))) for x in fault_traces[secid]])
    
    ''''''
    #---------------------------------------------------------------------------
    # t vs. P(t).
    #---------------------------------------------------------------------------
    cumulative = {}
    
    cumulative['x'] = np.sort(intervals)
    cumulative['y'] = np.arange(float(intervals.size))/float(intervals.size)
    #cumulative['y'] = [float(n)/float(len(intervals)) for n,x in enumerate(cumulative['x'])]
    #mplt.plot(np.sort(intervals), [float(n)/float(len(intervals)) for n,x in enumerate(np.sort(intervals))])
    
    #---------------------------------------------------------------------------
    # t0 vs. P(t0 + dt, t0) for fixed dt.
    #---------------------------------------------------------------------------
    
    conditional_dt_fixed = {'x':[],'y':[]}
    
    for t0 in np.sort(intervals):
        int_t0_dt = intervals[np.where( intervals > t0+fixed_dt)]
        int_t0 = intervals[np.where( intervals > t0)]
        
        if int_t0.size != 0:
            conditional_dt_fixed['x'].append(t0)
            conditional_dt_fixed['y'].append(1.0 - float(int_t0_dt.size)/float(int_t0.size))

    #---------------------------------------------------------------------------
    # t = t0 + dt vs. P(t, t0) for various t0.
    #---------------------------------------------------------------------------
    conditional = {}

    for t0 in range(0,275,25):
        int_t0 = intervals[np.where( intervals > t0)]
        if int_t0.size != 0:
            conditional[t0] = {'x':[],'y':[]}
            for dt in range(250):
                int_t0_dt = intervals[np.where( intervals > t0+dt)]
                conditional[t0]['x'].append(t0+dt)
                conditional[t0]['y'].append(1.0 - float(int_t0_dt.size)/float(int_t0.size))
        else:
            conditional[t0] = None

    #---------------------------------------------------------------------------
    # t0 vs dt.
    #---------------------------------------------------------------------------
    t0_dt = {}

    for percent in [0.25, 0.5, 0.75]:
        t0_dt[int(percent*100)] = {'x':[],'y':[]}
        for t0 in sorted(conditional.keys()):
            if conditional[t0] is not None:
                index = (np.abs(np.array(conditional[t0]['y'])-percent)).argmin()
                
                t0_dt[int(percent*100)]['x'].append(t0)
                t0_dt[int(percent*100)]['y'].append(conditional[t0]['x'][index]-t0)
        
    #---------------------------------------------------------------------------
    # Section probability matrix
    #---------------------------------------------------------------------------
    if event_graph_file is None:
        # Calculate the graph
        event_graph_file = 'event_graph.pkl'
        vcanalysis.graph_events(sim_file, event_graph_file, event_range=event_range, magnitude_filter=magnitude_filter)

    G = cPickle.load(open(event_graph_file, 'rb'))

    all_sections = [int(x) for x in sorted(fault_traces.keys())]

    matrix_prob = nx.attr_matrix(G, edge_attr='weight', normalized=False, rc_order=all_sections)
    
    if section_filter is not None:
        filtered_sections_pos = [np.where(np.array(all_sections) == sid)[0][0] for sid in sorted(section_filter['filter'])]
    
    #---------------------------------------------------------------------------
    # Sequences that end in the selected sections
    #---------------------------------------------------------------------------
    if event_sequence_graph_file is None:
        # Calculate the graph
        event_sequence_graph_file = 'event_sequence_graph.pkl'
        vcanalysis.graph_event_sequences(sim_file, event_sequence_graph_file, event_range=event_range, magnitude_filter=magnitude_filter)

    G = cPickle.load(open(event_sequence_graph_file, 'rb'))

    if section_filter is not None:
        sequences = {}
        for sid in sorted(section_filter['filter']):
            for sequence in G.to_undirected().neighbors_iter(sid):
                try:
                    sequences[sequence] += 1
                except KeyError:
                    sequences[sequence] = 1
    
    #---------------------------------------------------------------------------
    # Set up the plot.
    #---------------------------------------------------------------------------
    # The main map basemap instance
    mm = Basemap(
        llcrnrlon=min_lon,
        llcrnrlat=min_lat,
        urcrnrlon=max_lon,
        urcrnrlat=max_lat,
        lat_0=(max_lat+min_lat)/2.0,
        lon_0=(max_lon+min_lon)/2.0,
        resolution=map_res,
        projection=map_proj,
        suppress_ticks=True
    )
    
    #aspect = 0.974736393014
    aspect = mm.aspect
    
    # The main map height and width
    #mmw = imw - lm - rm - ftw - sphs
    mmh = math.ceil(mmw*aspect)
    
    # The sub map height and width
    smw = math.ceil((imw - lm - rm - 2.0*sphs)/3.0)
    smh = math.ceil(smw*aspect)
    
    # The total image height
    imh = tm + bm + mmh + 3.0*sph + 2.0*smh + 5.0*spvs

    # Create the figure instance
    imwi = imw/res
    imhi = imh/res
    fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
    
    #---------------------------------------------------------------------------
    # Plot the main map.
    #---------------------------------------------------------------------------
    mm.ax = fig.add_axes((lm/imw, (bm + 3.0*sph + 2.0*smh + 5.0*spvs)/imh, mmw/imw, mmh/imh))

    # Draw a frame around the map.
    mm.drawmapboundary(color=map_frame_color, linewidth=map_frame_width, fill_color=water_color)
    
    # Fill the continents and lakes.
    mm.fillcontinents(color=land_color, lake_color=water_color)
    
    # draw coastlines, edge of map.
    mm.drawcoastlines(color=coastline_color, linewidth=coastline_width)

    # draw countries
    mm.drawcountries(linewidth=country_width, color=country_color)

    # draw states
    mm.drawstates(linewidth=state_width, color=state_color)

    # draw parallels.
    parallels = np.linspace(min_lat, max_lat, num_grid_lines+1)
    mm_parallels = mm.drawparallels(
        parallels,
        labels=[1,0,0,0],
        color=grid_color,
        fontproperties=ticklabelfont,
        fmt='%.2f',
        linewidth=grid_width,
        dashes=[1, 10]
    )

    # draw meridians
    meridians = np.linspace(min_lon, max_lon, num_grid_lines+1)
    mm_meridians = mm.drawmeridians(
        meridians,
        labels=[0,0,1,0],
        color=grid_color,
        fontproperties=ticklabelfont,
        fmt='%.2f',
        linewidth=grid_width,
        dashes=[1, 10]
    )
    
    # print faults on lon-lat plot
    for sid, sec_trace in fault_traces_latlon.iteritems():
        trace_Xs, trace_Ys = mm(sec_trace[1], sec_trace[0])
        
        if section_filter is not None and sid in section_filter['filter']:
            linewidth = forecast_fault_width
        else:
            linewidth = fault_width

        mm.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=linewidth, solid_capstyle='round', solid_joinstyle='round')
        
    #---------------------------------------------------------------------------
    # Print the forecast title and info.
    #---------------------------------------------------------------------------
    ft_ax = fig.add_axes(((lm + mmw + sphs)/imw, (bm + 3.0*sph + 2.0*smh + 5.0*spvs)/imh, (imw - lm - rm - mmw - sphs)/imw, mmh/imh))

    ft_ax.text(0, mmh/mmh, 'Forecast sections', fontproperties=titlefont1, transform=ft_ax.transAxes, va='top')
    add_space = 0.0
    if section_filter['filter'] is not None:
        if len(section_filter['filter']) > 8:
            fs_text = ', '.join([str(x) for x in section_filter['filter'][0:9]])
            fs_text = fs_text + '\n' + ', '.join([str(x) for x in section_filter['filter'][9:-1]])
            add_space = 1.0
        else:
            fs_text = ', '.join([str(x) for x in section_filter['filter']])
    else:
        fs_text = 'all'
    ft_ax.text(0, (mmh - 1.0*forcast_title_line_spacing)/mmh, fs_text, fontproperties=titlefont2, transform=ft_ax.transAxes, va='top')
    
    ft_ax.text(0, (mmh - (3.0+add_space)*forcast_title_line_spacing)/mmh, 'Event magnitudes', fontproperties=titlefont1, transform=ft_ax.transAxes, va='top')
    if magnitude_filter is not None:
        em_text = magnitude_filter
    else:
        em_text = 'all'
    ft_ax.text(0, (mmh - (4.0+add_space)*forcast_title_line_spacing)/mmh, em_text, fontproperties=titlefont2, transform=ft_ax.transAxes, va='top')

    ft_ax.text(0, (mmh - (6.0+add_space)*forcast_title_line_spacing)/mmh, 'Simulation years', fontproperties=titlefont1, transform=ft_ax.transAxes, va='top')
    ft_ax.text(0, (mmh - (7.0+add_space)*forcast_title_line_spacing)/mmh, '{:0.0f}'.format(event_data['event_range_duration']), fontproperties=titlefont2, transform=ft_ax.transAxes, va='top')
    
    ft_ax.text(0, (mmh - (9.0+add_space)*forcast_title_line_spacing)/mmh, 'Number of events', fontproperties=titlefont1, transform=ft_ax.transAxes, va='top')
    ft_ax.text(0, (mmh - (10.0+add_space)*forcast_title_line_spacing)/mmh, '{}'.format(len(event_data['event_number'])), fontproperties=titlefont2, transform=ft_ax.transAxes, va='top')

    ft_ax.axis('off')
    
    ''''''
    #---------------------------------------------------------------------------
    # Plot t vs. P(t).
    #---------------------------------------------------------------------------
    cum_spw = math.ceil((imw - lm - rm - sphs)/2.0)
    cum_ax = fig.add_axes((lm/imw, (bm + 2.0*sph + 2.0*smh + 4.0*spvs)/imh, cum_spw/imw, sph/imh))

    cum_ax.plot(cumulative['x'], cumulative['y'], color=sp_line_color, linewidth=sp_line_width)
    
    # set the fonts for the tick labels
    for label in cum_ax.xaxis.get_ticklabels()+cum_ax.yaxis.get_ticklabels():
        label.set_fontproperties(ticklabelfont)
    
    cum_ax.set_ylim((0.0, 1.0))
    
    cum_ax.set_ylabel('P(t)', fontproperties=framelabelfont)
    cum_ax.set_xlabel('t [years]', fontproperties=framelabelfont)

    #---------------------------------------------------------------------------
    # Plot t0 vs. P(t0 + dt, t0) for fixed dt.
    #---------------------------------------------------------------------------
    cond_dt_fixed_spw = math.ceil((imw - lm - rm - sphs)/2.0)
    cond_dt_fixed_ax = fig.add_axes(((lm + cum_spw + sphs)/imw, (bm + 2.0*sph + 2.0*smh + 4.0*spvs)/imh, cum_spw/imw, sph/imh))

    cond_dt_fixed_ax.plot(conditional_dt_fixed['x'], conditional_dt_fixed['y'], color=sp_line_color, linewidth=sp_line_width)
    
    # set the fonts for the tick labels
    for label in cond_dt_fixed_ax.xaxis.get_ticklabels()+cond_dt_fixed_ax.yaxis.get_ticklabels():
        label.set_fontproperties(ticklabelfont)

    cond_dt_fixed_ax.set_ylim((0.0, 1.0))

    cond_dt_fixed_ax.set_ylabel(r'P(t$_0$ + {}, t$_0$)'.format(int(fixed_dt)), fontproperties=framelabelfont)
    cond_dt_fixed_ax.set_xlabel(r't$_0$ [years]', fontproperties=framelabelfont)

    #---------------------------------------------------------------------------
    # Plot t = t0 + dt vs. P(t, t0) for various t0.
    #---------------------------------------------------------------------------
    cond_spw = math.ceil((imw - lm - rm - sphs)/2.0)
    cond_ax = fig.add_axes((lm/imw, (bm + 1.0*sph + 2.0*smh + 3.0*spvs)/imh, cond_spw/imw, sph/imh))

    t0s_to_plot = [k for k in conditional.keys() if k <= 150]
    
    for t0 in sorted(t0s_to_plot):
        if conditional[t0] is not None:
            color = sp_line_colormap(float(t0)/float(max(t0s_to_plot)))
            cond_ax.plot(conditional[t0]['x'], conditional[t0]['y'], color=color, linewidth=sp_line_width, label='{}'.format(t0))

    # set the fonts for the tick labels
    for label in cond_ax.xaxis.get_ticklabels()+cond_ax.yaxis.get_ticklabels():
        label.set_fontproperties(ticklabelfont)

    cond_ax.set_xlim((0.0, max(conditional.keys())))
    cond_ax.set_ylim((0.0, 1.0))

    cond_ax.set_ylabel(r'P(t, t$_0$)', fontproperties=framelabelfont)
    cond_ax.set_xlabel(r't = t$_0$ + $\Delta$t [years]', fontproperties=framelabelfont)
    
    cond_ax.legend(title=r't$_0$=', prop=legendfont, loc=legend_loc)

    #---------------------------------------------------------------------------
    # Plot t0 vs dt.
    #---------------------------------------------------------------------------
    t0_dt_spw = math.ceil((imw - lm - rm - sphs)/2.0)
    t0_dt_ax = fig.add_axes(((lm + cond_spw + sphs)/imw, (bm + 1.0*sph + 2.0*smh + 3.0*spvs)/imh, t0_dt_spw/imw, sph/imh))

    percents = t0_dt.keys()
    t0_dt_ax.fill_between(t0_dt[min(percents)]['x'], t0_dt[min(percents)]['y'], y2=t0_dt[max(percents)]['y'], linewidth=0, facecolor=t0_dt_range_color)
    
    for percent in t0_dt.iterkeys():
        if percent == min(percents):
            linewidth = t0_dt_sub_line_width
            color = t0_dt_sub_line_color
            linestyle = '--'
        elif percent == max(percents):
            linewidth = t0_dt_sub_line_width
            color = t0_dt_sub_line_color
            linestyle = ':'
        else:
            linewidth = t0_dt_main_line_width
            color = t0_dt_main_line_color
            linestyle = '-'
        t0_dt_ax.plot(t0_dt[percent]['x'], t0_dt[percent]['y'], color=color, linewidth=linewidth, linestyle=linestyle, label='{}%'.format(percent))
    
    # set the fonts for the tick labels
    for label in t0_dt_ax.xaxis.get_ticklabels()+t0_dt_ax.yaxis.get_ticklabels():
        label.set_fontproperties(ticklabelfont)

    t0_dt_ax.set_ylabel(r'$\Delta$t [years]', fontproperties=framelabelfont)
    t0_dt_ax.set_xlabel(r't$_0$ [years]', fontproperties=framelabelfont)

    t0_dt_ax.legend(title='event prob.', prop=legendfont, loc=legend_loc, handlelength=5)

    #---------------------------------------------------------------------------
    # Plot the section probability matrix.
    #---------------------------------------------------------------------------
    prob_spw = imw - lm - rm
    prob_sph = sph - prob_cbh - prob_cbs
    prob_ax = fig.add_axes((lm/imw, (bm + 2.0*smh + 2.0*spvs + prob_cbh + prob_cbs)/imh, prob_spw/imw, prob_sph/imh))
    
    prob_ax.imshow(matrix_prob[filtered_sections_pos,:], aspect='auto', interpolation='none', cmap=matrix_cmap)

    prob_ax.set_xticks(range(0,len(all_sections),10))
    prob_ax.set_xticklabels([all_sections[x] for x in range(0,len(all_sections),10)])

    if section_filter is not None:
        prob_ax.set_yticks(range(len(section_filter['filter'])))
        prob_ax.set_yticklabels(section_filter['filter'])

    # set the fonts for the tick labels
    for label in prob_ax.xaxis.get_ticklabels()+prob_ax.yaxis.get_ticklabels():
        label.set_fontproperties(ticklabelfont)

    prob_ax.set_ylabel('to section', fontproperties=framelabelfont)
    prob_ax.set_xlabel('from section', fontproperties=framelabelfont)

    # Create the colorbar
    prob_cb_ax  = fig.add_axes((lm/imw, (bm + 2.0*smh + 2.0*spvs)/imh, prob_spw/imw, prob_cbh/imh))
    #norm = mcolor.Normalize(vmin=matrix_prob[filtered_sections_pos,:].min(), vmax=matrix_prob[filtered_sections_pos,:].max())
    cb = mcolorbar.ColorbarBase(prob_cb_ax, cmap=matrix_cmap,
               #norm=norm,
               orientation='horizontal')
    #ticks = map(float, [self.min_mag,(self.min_mag+ self.max_mag)/2.0, self.max_mag])
    cb.set_ticks((0,1))
    cb.set_ticklabels(('least likely','most likely'))
    
    # Style and cleanup the colorbar ticks
    for n, label in enumerate(prob_cb_ax.xaxis.get_ticklabels()):
        label.set_fontproperties(cbticklabelfont)
        if n == 0:
            label.set_ha('left')
        elif n == len(prob_cb_ax.xaxis.get_ticklabels()) - 1:
            label.set_ha('right')
    for line in prob_cb_ax.xaxis.get_ticklines():
        line.set_alpha(0)
    
    # Set the colorbar label
    #prob_cb_ax.set_xlabel('Probability',position=(0,0), ha='left', fontproperties=cbtitlefont)


    #---------------------------------------------------------------------------
    # Plot the event sequence maps.
    #---------------------------------------------------------------------------
    map_row = 2.0
    map_num = 0
    for sequence, val in sorted(sequences.iteritems(), key=itemgetter(1), reverse=True)[0:6]:
        
        sequence_list = [int(x) for x in sequence.split('->')]
        
        sm = Basemap(
            llcrnrlon=min_lon,
            llcrnrlat=min_lat,
            urcrnrlon=max_lon,
            urcrnrlat=max_lat,
            lat_0=(max_lat+min_lat)/2.0,
            lon_0=(max_lon+min_lon)/2.0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        
        if map_num%3 == 0:
            map_row -= 1.0
        
        sm.ax = fig.add_axes(((lm + map_num%3*(sphs+smw))/imw, (bm + map_row*(smh+spvs))/imh, smw/imw, smh/imh))

        # Draw a frame around the map.
        sm.drawmapboundary(color=map_frame_color, linewidth=map_frame_width, fill_color=water_color)
        
        # Fill the continents and lakes.
        sm.fillcontinents(color=seq_land_color, lake_color=water_color)
        
        # draw coastlines, edge of map.
        sm.drawcoastlines(color=coastline_color, linewidth=coastline_width-0.5)

        # draw countries
        sm.drawcountries(linewidth=country_width-0.5, color=country_color)

        # draw states
        sm.drawstates(linewidth=state_width-0.5, color=state_color)

        # print faults on lon-lat plot
        for sid, sec_trace in fault_traces_latlon.iteritems():
            trace_Xs, trace_Ys = sm(sec_trace[1], sec_trace[0])
            
            if section_filter is not None and sid in section_filter['filter']:
                linewidth = seq_fault_width_max
                color = sequence_cmap(1.0-float(5.0)/5.0)
            elif sid in sequence_list:
                seq_pos = sequence_list.index(sid)
                linewidth = vcutils.linear_interp(seq_pos, 0, 5, seq_fault_width_min, seq_fault_width_max)
                color = sequence_cmap(1.0-float(seq_pos)/5.0)
            else:
                linewidth = fault_width
                color = alt_fault_color

            sm.plot(trace_Xs, trace_Ys, color=color, linewidth=linewidth, solid_capstyle='round', solid_joinstyle='round')
        
        text_left = 0.0
        for n, sid in enumerate(sequence_list):
            if sid < 10:
                text_width = 25.0
            elif sid >=10 and sid<100:
                text_width = 30.0
            else:
                text_width = 35.0

            if n == 0 or n == 1:
                color = '#000000'
            else:
                color = '#ffffff'
            
            sm.ax.text(text_left/smw, -17.0/smh, r'{} $\rightarrow$ '.format(sid), ha='left', color=color, fontproperties=smtitlefont, bbox=dict(facecolor=sequence_cmap(1.0-float(n)/5.0), edgecolor='none'), transform=sm.ax.transAxes)
            text_left += text_width

        sm.ax.text(text_left/smw, -17.0/smh, r'{} $\ldots$ {}'.format(section_filter['filter'][0], section_filter['filter'][-1]), ha='left', color='#ffffff', fontproperties=smtitlefont, bbox=dict(facecolor=sequence_cmap(0), edgecolor='none'), transform=sm.ax.transAxes)
        
        #sm.ax.set_xlabel(r' $\rightarrow$ '.join(sequence.split('->')) + r' $\rightarrow$  $\times$' ,position=(0,0), ha='left', fontproperties=smtitlefont)

        map_num += 1
    
    '''
    #---------------------------------------------------------------------------
    # Print the section key.
    #---------------------------------------------------------------------------
    sk_ax = fig.add_axes(((lm + mmw + sphs)/imw, (bm + 3.0*sph + 2.0*smh + 5.0*spvs)/imh, (imw - lm - rm - mmw - sphs)/imw, mmh/imh))

    total_items = float(len(section_names))
    
    half = math.ceil(total_items/2.0)
    col = -0.5
    i = 0.0
    for sid, sname in section_names.iteritems():
        
        if i%half == 0:
            col += 0.5
            row = half
        #print i%half, col, row, half
        sk_ax.text(col, row/half, '{sid}: {sname}'.format(sid=sid,sname=sname), fontproperties=sectionkeyfont,transform=sk_ax.transAxes)
        row -= 1.0
        i += 1
    
    sk_ax.axis('off')
    '''
    if output_file is not None:
        # Get the plot format and save the file
        plot_format = output_file.split('.')[-1]
        if plot_format != 'png' and plot_format != 'pdf':
            raise vcexceptions.PlotFormatNotSupported(plot_format)
        else:
            fig.savefig(output_file, format=plot_format, dpi=res)
        #fig.savefig(output_file, format='png', dpi=res)


#-------------------------------------------------------------------------------
# Save each individual subplot from plot_forecast
#-------------------------------------------------------------------------------
def forecast_plots(sim_file, event_graph_file=None, event_sequence_graph_file=None, event_range=None, section_filter=None, magnitude_filter=None, padding=0.08, fixed_dt=30.0, fname_tag=None,weibull=False,
    tau=166.0,beta=1.5,year_eval=None,dt_and_t0=None):
    #---------------------------------------------------------------------------
    # Plot parameters.
    #---------------------------------------------------------------------------
    imw = 816.0 # The full image width.
    # Image height is set based on the sub plots.
    sph = 200.0 # The sub-plot height
    # Sub-plot width is based on the specific sub-plot.
    lm = 15.0
    rm = 15.0
    tm = 15.0
    bm = 20.0
    sphs = 40.0 # The sub-plot horizontal spacing.
    spvs = 40.0 # The sub-plot vertical spacing.
    
    mmw = 790.0
    
    res = 72.0
    map_res = 'i'
    map_proj = 'cyl'
    
    titlefont1 = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=12)
    titlefont2 = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=14, weight='bold')
    sectionkeyfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=10)
    ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=14)
    framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=14)
    mapticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=19)
    legendfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=12)
    smtitlefont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=12, weight='bold')
    cbticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=12)
    
    matrix_cmap = mplt.get_cmap('hot_r')
    sequence_cmap = mplt.get_cmap('autumn')
    
    water_color             = '#bed5ff'
    land_color              = '#ffffff'
    seq_land_color          = '#ffffff'
    boundary_color          = '#000000',
    coastline_color         = '#9a9a9a'
    country_color           = '#9a9a9a'
    state_color             = '#9a9a9a'
    fault_color             = '#000000'
    alt_fault_color         = '#737373'
    selected_fault_color    = '#FFFFFF'
    map_tick_color          = '#000000'
    map_frame_color         = '#000000'
    grid_color              = '#000000'
    cb_fontcolor            = '#000000'

    boundary_width          = 1.0
    coastline_width         = 1.0
    country_width           = 1.0
    state_width             = 1.0
    fault_width             = 0.5
    forecast_fault_width    = 6.0
    seq_fault_width_max     = 6.0
    seq_fault_width_min     = 3.0
    map_frame_width         = 1.0
    grid_width              = 0.5
    num_grid_lines          = 5
    
    sp_line_color           = '#000000'
    sp_line_colormap        = sequence_cmap
    sp_line_width           = 2.0
    
    t0_dt_main_line_color   = '#000000'
    t0_dt_sub_line_color    = '#737373'
    t0_dt_main_line_width   = 2.0
    t0_dt_sub_line_width    = 1.0
    t0_dt_range_color       = sequence_cmap(0.99)
    
    prob_cbh                = 20.0
    prob_cbs                = 40.0
    
    forcast_title_line_spacing = 20.0
    
    legend_loc='best'
    
    
    if section_filter is not None:
        if not isinstance( section_filter['filter'], ( int, long ) ):
            num_secs = len(section_filter['filter'])
        else:
            num_secs = 1
    else:
        sys.exit("***Error, must specify sections for a forecast.\n\t e.g. section_filter={'filter':(27,83)}\n")
     


    #---------------------------------------------------------------------------
    # Get the data.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)

        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        
        event_data = events.get_event_data(['event_number', 'event_year', 'event_magnitude', 'event_range_duration'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
    
        geometry = VCGeometry(sim_data)
        
        min_lat = geometry.min_lat
        max_lat = geometry.max_lat
        min_lon = geometry.min_lon
        max_lon = geometry.max_lon
        base_lat = geometry.base_lat
        base_lon = geometry.base_lon

        fault_traces = geometry.get_fault_traces()

        section_names = {sid:geometry.get_section_name(sid) for sid in fault_traces.iterkeys()}    



    intervals = np.array([   x - event_data['event_year'][n-1]
                    for n,x in enumerate(event_data['event_year'])
                    if n != 0
                ])
    
    avg_interval = np.average(intervals)
    max_interval = intervals.max()
    num_events   = len(event_data['event_year'])
    
    
    # Print out the relevant summary data
    sys.stdout.write('\n\n=======================================')
    sys.stdout.write('\n{} earthquakes \n'.format(num_events))
    sys.stdout.write('---------------------------------------\n')
    sys.stdout.write('avg interval\tmax interval\n')
    sys.stdout.write('{:.2f}\t\t{:.2f}\n'.format(avg_interval,max_interval))
    sys.stdout.write('---------------------------------------\n')

    
    # Calculate the lat-lon range based on the min-max and the padding
    lon_range = max_lon - min_lon
    lat_range = max_lat - min_lat
    max_range = max((lon_range, lat_range))
    min_lon = min_lon - lon_range*padding
    min_lat = min_lat - lat_range*padding
    max_lon = max_lon + lon_range*padding
    max_lat = max_lat + lat_range*padding

    
    
    # A conversion instance for doing the lat-lon to x-y conversions
    convert = quakelib.Conversion(base_lat, base_lon)
    
    # Convert the fault traces to lat-lon
    fault_traces_latlon = {}
    for secid in fault_traces.iterkeys():
         fault_traces_latlon[secid] = zip(*[(lambda y: (y.lat(),y.lon()))(convert.convert2LatLon(quakelib.Vec3(x[0], x[1], x[2]))) for x in fault_traces[secid]])
    
    ''''''
    # Don't want to plot recurrence times longer than those that actually occur
    #   and want to dynamically select the number of curves
    max_t0      = int(max_interval)
    
    # Want somewhere between 4 and 6 different t0 values for conditional plots
    if max_t0 < 100:
        t0_step = 5
        max_t0_cond = 20
        t0_plot_step = 3.0
    elif max_t0 >= 100 and max_t0 < 200:
        t0_step = 20
        max_t0_cond = 80
        t0_plot_step = 10.0
    elif max_t0 >= 200 and max_t0 < 250:
        t0_step = 30
        max_t0_cond = 120
        t0_plot_step = 20.0
    elif max_t0 >= 250 and max_t0 < 370:
        t0_step = 40
        max_t0_cond = 150
        t0_plot_step = 30.0
    else:
        t0_step      = 50
        max_t0_cond  = 250
        t0_plot_step = 40.0
        
    t0_for_multi_t0_plot = np.arange(0.0,max_t0_cond+t0_step*.1,t0_step)
    
    
    
    #---------------------------------------------------------------------------
    # t vs. P(t).
    #---------------------------------------------------------------------------
    cumulative = {}
    
    cumulative['x'] = np.sort(intervals)
    cumulative['y'] = np.arange(float(intervals.size))/float(intervals.size)
    
    #---------------------------------------------------------------------------
    # t0 vs. P(t0 + dt, t0) for fixed dt.
    #---------------------------------------------------------------------------
    
    conditional_dt_fixed = {'x':[],'y':[]}
    
    # To avoid the P(t0+30,t0) plot -> 0 at large t0, the probabilities are evaluated
    #  at incremental values up to the max interval, instead of evaluated at every
    #  interval
    t0_eval_step = 1.0
    t0_to_eval   = np.arange(0.0,max_t0+t0_eval_step*.01,t0_eval_step)

    
    for t0 in t0_to_eval:
        int_t0_dt = intervals[np.where( intervals > t0+fixed_dt)]
        int_t0 = intervals[np.where( intervals > t0)]
        
        if int_t0.size != 0:
            conditional_dt_fixed['x'].append(t0)
            conditional_dt_fixed['y'].append(1.0 - float(int_t0_dt.size)/float(int_t0.size))

    #---------------------------------------------------------------------------
    # t = t0 + dt vs. P(t, t0) for various t0.
    #---------------------------------------------------------------------------
    conditional = {}
    weibull     = {}

            
    #for t0 in range(0,max_t0,t0_step):
    for t0 in t0_to_eval:
        int_t0 = intervals[np.where( intervals > t0)]
        if int_t0.size != 0:
            conditional[t0] = {'x':[],'y':[]}
            weibull[t0]     = {'x':[],'y':[]}
            
            for dt in range(max_t0-int(t0)):
                int_t0_dt = intervals[np.where( intervals > t0+dt)]
                conditional[t0]['x'].append(t0+dt)
                weibull[t0]['x'].append(t0+dt)
                
                prob_t0_dt    = 1.0 - float(int_t0_dt.size)/float(int_t0.size)
                weibull_t0_dt = vcanalysis.cond_weibull(weibull[t0]['x'][-1],t0,beta,tau,single=True)   
                
                weibull[t0]['y'].append(weibull_t0_dt)
                conditional[t0]['y'].append(prob_t0_dt)
        else:
            conditional[t0] = None
            weibull[t0] = None
        


    #---------------------------------------------------------------------------
    # t0 vs dt.
    #---------------------------------------------------------------------------
    t0_dt      = {}
    t0_dt_plot = {}
    t0_dt_weibull = {}
    
    # In order to have a smooth plot, but still be able to evaluate the waiting
    #   times by specifying a year, I will have one t0_dt array for plotting
    #   (t0_dt_plot), and one for evaluating each year (t0_dt)
    

    for percent in [0.25, 0.5, 0.75]:
        t0_dt[int(percent*100)]      = {'x':[],'y':[]}
        t0_dt_plot[int(percent*100)] = {'x':[],'y':[]}
        
        t0_dt_weibull[int(percent*100)] = {'x':[],'y':[]}
        
        #for t0 in sorted(conditional.keys()):
        for t0 in t0_to_eval:
            if t0 < max_t0:
                if conditional[t0] is not None:
                    index   = (np.abs(np.array(conditional[t0]['y'])-percent)).argmin()
                    index_w = (np.abs(np.array(weibull[t0]['y'])-percent)).argmin()
                
                    t0_dt[int(percent*100)]['x'].append(t0)
                    t0_dt[int(percent*100)]['y'].append(conditional[t0]['x'][index]-t0)     
                     
                        
                    if int(t0)%int(t0_plot_step) == 0:
                        t0_dt_plot[int(percent*100)]['x'].append(t0)
                        t0_dt_plot[int(percent*100)]['y'].append(conditional[t0]['x'][index]-t0)
                        
                        t0_dt_weibull[int(percent*100)]['x'].append(t0)
                        t0_dt_weibull[int(percent*100)]['y'].append(weibull[t0]['x'][index_w]-t0)
                
    


    #---------------------------------------------------------------------------
    # Grab specific recurrence times
    #---------------------------------------------------------------------------
    if year_eval is not None:

        wait_25  = None
        wait_50  = None
        wait_75  = None
        
        ind_25 = (np.abs(np.array(t0_dt[25]['x'])-year_eval)).argmin()
        ind_50 = (np.abs(np.array(t0_dt[50]['x'])-year_eval)).argmin()
        ind_75 = (np.abs(np.array(t0_dt[75]['x'])-year_eval)).argmin()

        wait_25 = t0_dt[25]['y'][ind_25]        
        wait_50 = t0_dt[50]['y'][ind_50]
        wait_75 = t0_dt[75]['y'][ind_75]
                    
        sys.stdout.write('For t0 = {:.2f} years'.format(year_eval))
        sys.stdout.write('\n25% waiting time: {:.2f} years'.format(wait_25))
        sys.stdout.write('\n50% waiting time: {:.2f} years'.format(wait_50))
        sys.stdout.write('\n75% waiting time: {:.2f} years'.format(wait_75))
        sys.stdout.write('\n=======================================\n\n')
    
    if dt_and_t0 is not None:
        dt = int(dt_and_t0[0])
        t0 = int(dt_and_t0[1])

        # P(t0 < t < t0+dt) = cond[t0+dt] - cond[t0]
        #upper = (np.abs(np.array(conditional[t0]['x'])-t0+dt)).argmin()
        #lower = (np.abs(np.array(conditional[t0]['x'])-t0)).argmin()
        #prob_in_t0_dt = conditional[t0]['y'][upper] - conditional[t0]['y'][lower]

        prob_in_t0_dt = conditional[t0]['y'][dt] - conditional[t0]['y'][0]

        #sys.stdout.write('\n upper/lower index for t0/t0+dt={}/{} is {}/{} \n'.format(t0,t0+dt,lower,upper))

        sys.stdout.write('\n ---------- evaluating P(t0 < t < t0+dt) -------')
        sys.stdout.write('\n P({} < t < {}) = {}\n\n'.format(t0,t0+dt,prob_in_t0_dt))



        
    """    
    #---------------------------------------------------------------------------
    # Section probability matrix
    #---------------------------------------------------------------------------
    if event_graph_file is None:
        # Calculate the graph
        event_graph_file = 'event_graph.pkl'
        vcanalysis.graph_events(sim_file, event_graph_file, event_range=event_range, magnitude_filter=magnitude_filter)

    G = cPickle.load(open(event_graph_file, 'rb'))

    all_sections = [int(x) for x in sorted(fault_traces.keys())]

    matrix_prob = nx.attr_matrix(G, edge_attr='weight', normalized=False, rc_order=all_sections)
    
    if section_filter is not None:
        filtered_sections_pos = [np.where(np.array(all_sections) == sid)[0][0] for sid in sorted(section_filter['filter'])]
    
    #---------------------------------------------------------------------------
    # Sequences that end in the selected sections
    #---------------------------------------------------------------------------
    if event_sequence_graph_file is None:
        # Calculate the graph
        event_sequence_graph_file = 'event_sequence_graph.pkl'
        vcanalysis.graph_event_sequences(sim_file, event_sequence_graph_file, event_range=event_range, magnitude_filter=magnitude_filter)

    G = cPickle.load(open(event_sequence_graph_file, 'rb'))

    if section_filter is not None:
        sequences = {}
        for sid in sorted(section_filter['filter']):
            for sequence in G.to_undirected().neighbors_iter(sid):
                try:
                    sequences[sequence] += 1
                except KeyError:
                    sequences[sequence] = 1
    """
    
    #---------------------------------------------------------------------------
    # Set up the plot.
    #---------------------------------------------------------------------------
    # The main map basemap instance
    mm = Basemap(
        llcrnrlon=min_lon,
        llcrnrlat=min_lat,
        urcrnrlon=max_lon,
        urcrnrlat=max_lat,
        lat_0=(max_lat+min_lat)/2.0,
        lon_0=(max_lon+min_lon)/2.0,
        resolution=map_res,
        projection=map_proj,
        suppress_ticks=True
    )
    
    #aspect = 0.974736393014
    aspect = mm.aspect
    
    # The main map height and width
    #mmw = imw - lm - rm - ftw - sphs
    mmh = math.ceil(mmw*aspect)
    
    # The sub map height and width
    smw = math.ceil((imw - lm - rm - 2.0*sphs)/3.0)
    smh = math.ceil(smw*aspect)
    
    # The total image height
    imh = tm + bm + mmh + 3.0*sph + 2.0*smh + 5.0*spvs

    #---------------------------------------------------------------------------
    # Create the figure instance for big map
    #---------------------------------------------------------------------------
    imwi = mmw/res
    imhi = mmh/res
    fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
    
    #---------------------------------------------------------------------------
    # Plot the main map.
    #---------------------------------------------------------------------------
    #mm.ax = fig.add_axes((lm/imw, (bm + 3.0*sph + 2.0*smh + 5.0*spvs)/imh, mmw/imw, mmh/imh))
    mm.ax = fig.add_axes()

    # Draw a frame around the map.
    mm.drawmapboundary(color=map_frame_color, linewidth=map_frame_width, fill_color=water_color)
    
    # Fill the continents and lakes.
    mm.fillcontinents(color=land_color, lake_color=water_color)
    
    # draw coastlines, edge of map.
    mm.drawcoastlines(color=coastline_color, linewidth=coastline_width)

    # draw countries
    mm.drawcountries(linewidth=country_width, color=country_color)

    # draw states
    mm.drawstates(linewidth=state_width, color=state_color)

    # draw parallels.
    parallels = np.linspace(min_lat, max_lat, num_grid_lines+1)
    mm_parallels = mm.drawparallels(
        parallels,
        labels=[1,0,0,0],
        color=grid_color,
        fontproperties=mapticklabelfont,
        fmt='%.2f',
        linewidth=grid_width,
        dashes=[1, 10]
    )

    # draw meridians
    meridians = np.linspace(min_lon, max_lon, num_grid_lines+1)
    mm_meridians = mm.drawmeridians(
        meridians,
        labels=[0,0,1,0],
        color=grid_color,
        fontproperties=mapticklabelfont,
        fmt='%.2f',
        linewidth=grid_width,
        dashes=[1, 10]
    )
    
    # print faults on lon-lat plot
    for sid, sec_trace in fault_traces_latlon.iteritems():
        trace_Xs, trace_Ys = mm(sec_trace[1], sec_trace[0])
        
        if section_filter is not None and sid in section_filter['filter']:
            linewidth = forecast_fault_width
        else:
            linewidth = fault_width

        mm.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=linewidth, solid_capstyle='round', solid_joinstyle='round')
    
    fname1 = 'forecast_sections.png'
    
    if fname_tag is not None:
        outfile1 = 'local/'+fname_tag+'_'+fname1
    else:
        outfile1 = 'local/'+fname1
                
    mplt.savefig(outfile1, dpi=res*3)
        
    
    ''''''
    #---------------------------------------------------------------------------
    # Plot t vs. P(t).
    #---------------------------------------------------------------------------

    imwi = imw/(res*1.5)
    imhi = sph*2.0/res
    
    fig2    = mplt.figure(figsize=(imwi, imhi), dpi=res)
    cum_ax  = fig2.add_subplot(111)
    
    #  Plot the weibull distribution -------------------------------------------
    if weibull:
        cum_ax.plot(cumulative['x'], vcanalysis.weibull(cumulative['x'],beta,tau),color='r',linewidth=sp_line_width-1.0)


    cum_ax.plot(cumulative['x'], cumulative['y'], color=sp_line_color, linewidth=sp_line_width)
    
    # set the fonts for the tick labels
    for label in cum_ax.xaxis.get_ticklabels()+cum_ax.yaxis.get_ticklabels():
        label.set_fontproperties(ticklabelfont)
    
    cum_ax.set_ylim((0.0, 1.0))
    
    cum_ax.set_ylabel('P(t)', fontproperties=framelabelfont)
    cum_ax.set_xlabel('t [years]', fontproperties=framelabelfont)
    
    fname2 = 'prob_vs_time.png'
    
    if fname_tag is not None:
        outfile2 = 'local/'+fname_tag+'_'+fname2
    else:
        outfile2 = 'local/'+fname2
                
    mplt.savefig(outfile2, dpi=res*2)

    #---------------------------------------------------------------------------
    # Plot t0 vs. P(t0 + dt, t0) for fixed dt.
    #---------------------------------------------------------------------------
    fig3 = mplt.figure(figsize=(imwi, imhi), dpi=res)
    
    #cond_dt_fixed_spw = math.ceil((imw - lm - rm - sphs)/2.0)
    #cond_dt_fixed_ax = fig.add_axes(((lm + cum_spw + sphs)/imw, (bm + 2.0*sph + 2.0*smh + 4.0*spvs)/imh, cum_spw/imw, sph/imh))
    cond_dt_fixed_ax = fig3.add_subplot(111)

    if weibull:
        weibull_vals = vcanalysis.cond_weibull_fixed_dt(conditional_dt_fixed['x'],fixed_dt,beta,tau)
        cond_dt_fixed_ax.plot(conditional_dt_fixed['x'],weibull_vals, color='r', linewidth=sp_line_width-1.0)

    cond_dt_fixed_ax.plot(conditional_dt_fixed['x'], conditional_dt_fixed['y'], color=sp_line_color, linewidth=sp_line_width)
    
    
    # set the fonts for the tick labels
    for label in cond_dt_fixed_ax.xaxis.get_ticklabels()+cond_dt_fixed_ax.yaxis.get_ticklabels():
        label.set_fontproperties(ticklabelfont)

    cond_dt_fixed_ax.set_ylim((0.0, 1.0))

    cond_dt_fixed_ax.set_ylabel(r'P(t$_0$ + {}, t$_0$)'.format(int(fixed_dt)), fontproperties=framelabelfont)
    cond_dt_fixed_ax.set_xlabel(r't$_0$ [years]', fontproperties=framelabelfont)
    
    fname3 = 'cond_prob_'+str(int(fixed_dt))+'_vs_time.png'
    
    if fname_tag is not None:
        outfile3 = 'local/'+fname_tag+'_'+fname3
    else:
        outfile3 = 'local/'+fname3
                
    mplt.savefig(outfile3, dpi=res*2)


    #---------------------------------------------------------------------------
    # Plot t = t0 + dt vs. P(t, t0) for various t0.
    #---------------------------------------------------------------------------
    fig4 = mplt.figure(figsize=(imwi, imhi), dpi=res)
    
    cond_ax = fig4.add_subplot(111)

    #t0s_to_plot = [k for k in conditional.keys() if k <= max_t0_cond]
    t0s_to_plot = t0_for_multi_t0_plot
    
    for t0 in sorted(t0s_to_plot):
        if conditional[t0] is not None:
            color = sp_line_colormap(float(t0*.8)/float(max(t0s_to_plot)))

            if weibull:
                cond_ax.plot(conditional[t0]['x'], vcanalysis.cond_weibull(conditional[t0]['x'],t0,beta,tau), color='k', linewidth=sp_line_width-1.0)

            cond_ax.plot(conditional[t0]['x'], conditional[t0]['y'], color=color, linewidth=sp_line_width, label='{}'.format(int(t0)))



    # set the fonts for the tick labels
    for label in cond_ax.xaxis.get_ticklabels()+cond_ax.yaxis.get_ticklabels():
        label.set_fontproperties(ticklabelfont)

    cond_ax.set_xlim((0.0, max(conditional.keys())))
    cond_ax.set_ylim((0.0, 1.0))

    cond_ax.set_ylabel(r'P(t, t$_0$)', fontproperties=framelabelfont)
    cond_ax.set_xlabel(r't = t$_0$ + $\Delta$t [years]', fontproperties=framelabelfont)
    
    cond_ax.legend(title=r't$_0$=', prop=legendfont, loc=legend_loc)
    
    fname4 = 'cum_prob_vs_multi_t0.png'
    
    if fname_tag is not None:
        outfile4 = 'local/'+fname_tag+'_'+fname4
    else:
        outfile4 = 'local/'+fname4
                
    mplt.savefig(outfile4, dpi=res*2)


    #---------------------------------------------------------------------------
    # Plot t0 vs dt.
    #---------------------------------------------------------------------------
    fig5 = mplt.figure(figsize=(imwi, imhi), dpi=res)
    
    #t0_dt_spw = math.ceil((imw - lm - rm - sphs)/2.0)
    #t0_dt_ax = fig.add_axes(((lm + cond_spw + sphs)/imw, (bm + 1.0*sph + 2.0*smh + 3.0*spvs)/imh, t0_dt_spw/imw, sph/imh))
    t0_dt_ax = fig5.add_subplot(111)

    percents = t0_dt.keys()
    t0_dt_ax.fill_between(t0_dt_plot[min(percents)]['x'], t0_dt_plot[min(percents)]['y'], y2=t0_dt_plot[max(percents)]['y'], linewidth=0, facecolor=t0_dt_range_color)
    
    for percent in t0_dt_plot.iterkeys():
        if percent == min(percents):
            linewidth = t0_dt_sub_line_width
            color = t0_dt_sub_line_color
            linestyle = '--'
        elif percent == max(percents):
            linewidth = t0_dt_sub_line_width
            color = t0_dt_sub_line_color
            linestyle = ':'
        else:
            linewidth = t0_dt_main_line_width
            color = t0_dt_main_line_color
            linestyle = '-'

        #if weibull:
        #    t0_dt_ax.plot(t0_dt_weibull[percent]['x'], t0_dt_weibull[percent]['y'], color='r', linewidth=t0_dt_sub_line_width, linestyle='-')


        t0_dt_ax.plot(t0_dt_plot[percent]['x'], t0_dt_plot[percent]['y'], color=color, linewidth=linewidth, linestyle=linestyle, label='{}%'.format(percent))
        

        
    
    if year_eval is not None and year_eval != 0:
        t0_dt_ax.axvline(x=year_eval,ymin=0,ymax=wait_75,color='blue',linewidth=t0_dt_main_line_width,linestyle='--')
    
    
    # set the fonts for the tick labels
    for label in t0_dt_ax.xaxis.get_ticklabels()+t0_dt_ax.yaxis.get_ticklabels():
        label.set_fontproperties(ticklabelfont)

    t0_dt_ax.set_ylabel(r'$\Delta$t [years]', fontproperties=framelabelfont)
    t0_dt_ax.set_xlabel(r't$_0$ [years]', fontproperties=framelabelfont)

    t0_dt_ax.legend(title='event prob.', prop=legendfont, loc=legend_loc, handlelength=5)
    
    
    fname5 = 't0_vs_dt.png'
    
    if fname_tag is not None:
        outfile5 = 'local/'+fname_tag+'_'+fname5
    else:
        outfile5 = 'local/'+fname5
                
    mplt.savefig(outfile5, dpi=res*2)

    
    """
    #---------------------------------------------------------------------------
    # Plot the section probability matrix.
    #---------------------------------------------------------------------------
    #prob_spw = imw - lm - rm
    #prob_sph = sph - prob_cbh - prob_cbs
    # More forecast sections => more height needed in matrix
    mat_height = (1.0 + num_secs/5)*80.0
    
    imwi    = imw/res
    imhi    = mat_height/res
    fig6    = mplt.figure(figsize=(imwi, imhi), dpi=res)
    prob_ax = fig6.add_axes((0.075,.21,.9,.78))
    
    #prob_ax = fig.add_axes((lm/imw, (bm + 2.0*smh + 2.0*spvs + prob_cbh + prob_cbs)/imh, prob_spw/imw, prob_sph/imh))
    
    prob_ax.imshow(matrix_prob[filtered_sections_pos,:], aspect='auto', interpolation='none', cmap=matrix_cmap)

    prob_ax.set_xticks(range(0,len(all_sections),10))
    prob_ax.set_xticklabels([all_sections[x] for x in range(0,len(all_sections),10)])

    if section_filter is not None:
        prob_ax.set_yticks(range(len(section_filter['filter'])))
        prob_ax.set_yticklabels(section_filter['filter'])

    # set the fonts for the tick labels
    for label in prob_ax.xaxis.get_ticklabels()+prob_ax.yaxis.get_ticklabels():
        label.set_fontproperties(ticklabelfont)

    prob_ax.set_ylabel('to section', fontproperties=framelabelfont)
    prob_ax.set_xlabel('from section', fontproperties=framelabelfont)

    # Create the colorbar
    #prob_cb_ax  = fig.add_axes((lm/imw, (bm + 2.0*smh + 2.0*spvs)/imh, prob_spw/imw, prob_cbh/imh))
    prob_cb_ax   = fig6.add_axes((0.075,0.05,.9,.05))
    #norm = mcolor.Normalize(vmin=matrix_prob[filtered_sections_pos,:].min(), vmax=matrix_prob[filtered_sections_pos,:].max())
    cb = mcolorbar.ColorbarBase(prob_cb_ax, cmap=matrix_cmap,
               #norm=norm,
               orientation='horizontal')
    #ticks = map(float, [self.min_mag,(self.min_mag+ self.max_mag)/2.0, self.max_mag])
    cb.set_ticks((0,1))
    cb.set_ticklabels(('least likely','most likely'))
    
    # Style and cleanup the colorbar ticks
    for n, label in enumerate(prob_cb_ax.xaxis.get_ticklabels()):
        label.set_fontproperties(cbticklabelfont)
        if n == 0:
            label.set_ha('left')
        elif n == len(prob_cb_ax.xaxis.get_ticklabels()) - 1:
            label.set_ha('right')
    for line in prob_cb_ax.xaxis.get_ticklines():
        line.set_alpha(0)
    
    # Set the colorbar label
    #prob_cb_ax.set_xlabel('Probability',position=(0,0), ha='left', fontproperties=cbtitlefont)

    fname6 = 'trans_prob_matrix.png'
    
    if fname_tag is not None:
        outfile6 = 'local/'+fname_tag+'_'+fname6
    else:
        outfile6 = 'local/'+fname6
                
    mplt.savefig(outfile6, dpi=res)

    
    #---------------------------------------------------------------------------
    # Plot the event sequence maps.
    #---------------------------------------------------------------------------
    imw     = 816.0
    imh     = 544.0
    imwi    = imw/res
    imhi    = imh/res
    fig6    = mplt.figure(figsize=(imwi, imhi), dpi=res)

    map_row = 2.0
    map_num = 0
    for sequence, val in sorted(sequences.iteritems(), key=itemgetter(1), reverse=True)[0:6]:
        
        sequence_list = [int(x) for x in sequence.split('->')]
        
        sm = Basemap(
            llcrnrlon=min_lon,
            llcrnrlat=min_lat,
            urcrnrlon=max_lon,
            urcrnrlat=max_lat,
            lat_0=(max_lat+min_lat)/2.0,
            lon_0=(max_lon+min_lon)/2.0,
            resolution=map_res,
            projection=map_proj,
            suppress_ticks=True
        )
        
        if map_num%3 == 0:
            map_row -= 1.0
        
        left   = (lm + map_num%3*(sphs+smw))/imw
        bottom = (bm + map_row*(smh+spvs))/imh
        width  = smw/imw
        height = smh/imh
        sm.ax  = fig6.add_axes((left, bottom, width, height))

        # Draw a frame around the map.
        sm.drawmapboundary(color=map_frame_color, linewidth=map_frame_width, fill_color=water_color)
        
        # Fill the continents and lakes.
        sm.fillcontinents(color=seq_land_color, lake_color=water_color)
        
        # draw coastlines, edge of map.
        sm.drawcoastlines(color=coastline_color, linewidth=coastline_width-0.5)

        # draw countries
        sm.drawcountries(linewidth=country_width-0.5, color=country_color)

        # draw states
        sm.drawstates(linewidth=state_width-0.5, color=state_color)

        # print faults on lon-lat plot
        for sid, sec_trace in fault_traces_latlon.iteritems():
            trace_Xs, trace_Ys = sm(sec_trace[1], sec_trace[0])
            
            if section_filter is not None and sid in section_filter['filter']:
                linewidth = seq_fault_width_max
                color = sequence_cmap(1.0-float(5.0)/5.0)
            elif sid in sequence_list:
                seq_pos = sequence_list.index(sid)
                linewidth = vcutils.linear_interp(seq_pos, 0, 5, seq_fault_width_min, seq_fault_width_max)
                color = sequence_cmap(1.0-float(seq_pos)/5.0)
            else:
                linewidth = fault_width
                color = alt_fault_color

            sm.plot(trace_Xs, trace_Ys, color=color, linewidth=linewidth, solid_capstyle='round', solid_joinstyle='round')
        
        text_left = 0.0
        for n, sid in enumerate(sequence_list):
            if sid < 10:
                text_width = 25.0
            elif sid >=10 and sid<100:
                text_width = 30.0
            else:
                text_width = 35.0

            if n == 0 or n == 1:
                color = '#000000'
            else:
                color = '#ffffff'
            
            sm.ax.text(text_left/smw, -17.0/smh, r'{} $\rightarrow$ '.format(sid), ha='left', color=color, fontproperties=smtitlefont, bbox=dict(facecolor=sequence_cmap(1.0-float(n)/5.0), edgecolor='none'), transform=sm.ax.transAxes)
            text_left += text_width

        sm.ax.text(text_left/smw, -17.0/smh, r'{} $\ldots$ {}'.format(section_filter['filter'][0], section_filter['filter'][-1]), ha='left', color='#ffffff', fontproperties=smtitlefont, bbox=dict(facecolor=sequence_cmap(0), edgecolor='none'), transform=sm.ax.transAxes)
        
        #sm.ax.set_xlabel(r' $\rightarrow$ '.join(sequence.split('->')) + r' $\rightarrow$  $\times$' ,position=(0,0), ha='left', fontproperties=smtitlefont)

        map_num += 1
    
    
    
    
    fname7 = 'hi_prob_sequences.png'
    
    if fname_tag is not None:
        outfile7 = 'local/'+fname_tag+'_'+fname7
    else:
        outfile7 = 'local/'+fname7
                
    mplt.savefig(outfile7, dpi=res)
    """




#-------------------------------------------------------------------------------
# event field animation
#-------------------------------------------------------------------------------
def event_field_animation(sim_file, output_directory, event_range,
    field_type='displacement', fringes=True, padding=0.08, cutoff=None,
    animation_target_length=60.0, animation_fps = 30.0, fade_seconds = 1.0,
    min_mag_marker = 6.5, force_plot=False, contours=None):
    
    sys.stdout.write('Initializing animation :: ')
    sys.stdout.flush()
    
    # create the animation dir if needed
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    if not output_directory.endswith('/'):
        output_directory += '/'
    
    # create the animation subdirs if needed
    field_values_directory = '{}field_values/'.format(output_directory)
    frame_images_directory = '{}frame_images/'.format(output_directory)
    if not os.path.exists(field_values_directory):
        os.makedirs(field_values_directory)
    if not os.path.exists(frame_images_directory):
        os.makedirs(frame_images_directory)

    # animation properties

    #---------------------------------------------------------------------------
    # Open the data file. It needs to stay open while we do the animation.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)
        
        # Get global information about the simulations geometry
        min_lat = geometry.min_lat
        max_lat = geometry.max_lat
        min_lon = geometry.min_lon
        max_lon = geometry.max_lon
        base_lat = geometry.base_lat
        base_lon = geometry.base_lon
        fault_traces = geometry.get_fault_traces()
        
        # Get event information
        event_data = events.get_event_data(['event_magnitude', 'event_year', 'event_number'], event_range=event_range)
        event_magnitudes = event_data['event_magnitude']
        event_years = event_data['event_year']
        event_numbers = event_data['event_number']
        current_year = start_year = math.floor(event_years[0])
        
        # These are the event years shifted so the begin at zero.
        _event_years = [y - start_year for y in event_years]
        
        # The large magnitudes and shifted years to be marked on the timeline.
        event_large_magnitudes = [m for m in event_magnitudes if m > min_mag_marker]
        _event_large_magnitude_years = [_event_years[i_m[0]] for i_m in enumerate(event_magnitudes) if i_m[1] >= min_mag_marker]
        event_large_magnitude_evnums = [event_numbers[i_m[0]] for i_m in enumerate(event_magnitudes) if i_m[1] > min_mag_marker]
        
        # Calculate the frames per year and the total number of frames
        total_years = math.ceil(event_years[-1]) - math.floor(event_years[0])
        fpy = math.ceil(animation_target_length*animation_fps/total_years)
        total_frames = int(fpy * total_years)
        
        # Instantiate the field and the plotter
        if field_type == 'displacement':
            EF = vcutils.VCDisplacementField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
            EFP = vcplotutils.VCDisplacementFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)
            EFP.calculate_look_angles(geometry[:])
        elif field_type == 'gravity':
            EF = vcutils.VCGravityField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding,  free_air=True)
            EFP = vcplotutils.VCGravityFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon, contours=contours)
        elif field_type == 'dilat_gravity':
            EF = vcutils.VCGravityField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding,  free_air=False)
            EFP = vcplotutils.VCGravityFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon, contours=contours)

        #-----------------------------------------------------------------------
        # Find the biggest event and normalize based on these values.
        #-----------------------------------------------------------------------
        if field_type == 'displacement' and not fringes or field_type == 'gravity' or field_type == 'dilat_gravity':
            sys.stdout.write('normalizing : ')
            sys.stdout.flush()
            max_mag_evnum = event_numbers[event_magnitudes.index(max(event_magnitudes))]
            field_values_loaded = EF.load_field_values('{}{}_'.format(field_values_directory, max_mag_evnum))
            if field_values_loaded:
                sys.stdout.write('event {} loaded : '.format(max_mag_evnum))
                sys.stdout.flush()
            if not field_values_loaded:
                sys.stdout.write('event {} processing : '.format(max_mag_evnum))
                sys.stdout.flush()
                event_element_slips = events.get_event_element_slips(max_mag_evnum)
                ele_getter = itemgetter(*event_element_slips.keys())
                event_element_data = ele_getter(geometry)
                if len(event_element_slips) == 1:
                    event_element_data = [event_element_data]
            
                sys.stdout.write('{} elements : '.format(len(event_element_slips)))
                sys.stdout.flush()
                
                EF.calculate_field_values(
                    event_element_data,
                    event_element_slips,
                    cutoff=cutoff,
                    save_file_prefix='{}{}_'.format(field_values_directory,max_mag_evnum)
                )
            EFP.set_field(EF)
            EFP.create_field_image()
        
        # Convert the fault traces to lat-lon
        fault_traces_latlon = {}
        for secid in fault_traces.iterkeys():
             fault_traces_latlon[secid] = zip(*[(lambda y: (y.lat(),y.lon()))(EF.convert.convert2LatLon(quakelib.Vec3(x[0], x[1], x[2]))) for x in fault_traces[secid]])

        # Grab all of the plot properties that we will need.
        # properties that are fringes dependent
        if fringes and field_type == 'displacement':
            cmap            = EFP.dmc['cmap_f']
            coastline_color = EFP.dmc['coastline_color_f']
            country_color   = EFP.dmc['country_color_f']
            state_color     = EFP.dmc['state_color_f']
            fault_color     = EFP.dmc['fault_color_f']
            map_tick_color  = EFP.dmc['map_tick_color_f']
            map_frame_color = EFP.dmc['map_frame_color_f']
            grid_color      = EFP.dmc['grid_color_f']
            cb_fontcolor    = EFP.dmc['cb_fontcolor_f']
        else:
            cmap            = EFP.dmc['cmap']
            coastline_color = EFP.dmc['coastline_color']
            country_color   = EFP.dmc['country_color']
            state_color     = EFP.dmc['state_color']
            fault_color     = EFP.dmc['fault_color']
            map_tick_color  = EFP.dmc['map_tick_color']
            map_frame_color = EFP.dmc['map_frame_color']
            grid_color      = EFP.dmc['grid_color']
            cb_fontcolor    = EFP.dmc['cb_fontcolor']
        
        # properties that are not fringes dependent
        boundary_width  = EFP.dmc['boundary_width']
        coastline_width = EFP.dmc['coastline_width']
        country_width   = EFP.dmc['country_width']
        state_width     = EFP.dmc['state_width']
        river_width     = EFP.dmc['river_width']
        fault_width     = EFP.dmc['fault_width']
        map_frame_width = EFP.dmc['map_frame_width']
        map_fontsize    = EFP.dmc['map_fontsize']
        arrow_inset     = EFP.dmc['arrow_inset']
        arrow_fontsize  = EFP.dmc['arrow_fontsize']
        cb_fontsize     = EFP.dmc['cb_fontsize']
        cb_height       = EFP.dmc['cb_height']
        cb_margin_t     = EFP.dmc['cb_margin_t']
        grid_width      = EFP.dmc['grid_width']
        num_grid_lines  = EFP.dmc['num_grid_lines']
        font            = EFP.dmc['font']
        font_bold       = EFP.dmc['font_bold']

        map_resolution  = EFP.dmc['map_resolution']
        map_projection  = EFP.dmc['map_projection']
        plot_resolution = EFP.dmc['plot_resolution']

        #animation specific properties
        progress_tick_color = 'k'
        progress_frame_color = 'k'
        progress_frame_width = 1
        progress_line_color = 'k'
        progress_line_width = 0
        progress_marker = '.'
        progress_marker_edge_width = 0
        progress_marker_size = 4
        progress_indicator_line_color = 'red'
        progress_indicator_linewidth = 2
        progress_indicator_fontsize = 10.0
        
        mag_color = 'k'
        current_mag_color = 'white'
        current_mag_facecolor = 'red'
        mag_facecolor = 'white'
        mag_linewidth = 0.5
        mag_fontsize = 12
        
        fm_fontsize = 9
        if field_type=='displacement':
            fm_label_color = 'white'
        else:
            fm_label_color = 'black'
            
        fm_frame_width = 1
        fm_frame_color = 'k'
        fm_line_color = '0.0'
        fm_line_width = 1
        
        mag_color_map = mcolor.LinearSegmentedColormap.from_list(
            'mag_color_map',
            [mag_color,current_mag_color],
            N=256,
            gamma=1.0
        )
        mag_line_colormap = mcolor.LinearSegmentedColormap.from_list(
            'mag_line_colormap',
            [progress_frame_color,progress_indicator_line_color],
            N=256,
            gamma=1.0
        )
        current_mag_face_colormap = mcolor.LinearSegmentedColormap.from_list(
            'current_mag_face_colormap',
            [mag_facecolor,current_mag_facecolor],
            N=256,
            gamma=1.0
        )

        # Set up the field values. In order to do smooth fading, we are gonna be
        # adding to these quantities when there are new events and subtracting a
        # bit each frame till we get back down to zero.
        EF.init_field(0.0)

        # We need to keep track of all of the magnitudes that have been plotted
        # so far for the FM plot. Also, set up some other properties of the FM
        # plot.
        cumulative_magnitudes = []
        fm_x_ticks = np.linspace(round(min(event_magnitudes)), round(max(event_magnitudes)), 5)
        fm_y_ticks = np.logspace(0, 3, 4)
        
        # Set up a cuple of dicts to maintain the state of the large magnitude
        # labels and the section lines. This will allow us to fade these items
        # out after an event.
        section_states = {}
        large_event_label_states = {}
        fm_alpha_state = 0.0
        
        sys.stdout.write('done\n')
        sys.stdout.flush()
        
        #-----------------------------------------------------------------------
        # Go through all of the frames.
        #-----------------------------------------------------------------------
        sys.stdout.write('Total frames : {}, Frames per year : {}\n'.format(total_frames, fpy))
        sys.stdout.flush()
        for the_frame in range(total_frames):
            year_frame = the_frame%fpy
            if year_frame == 0:
                current_year += 1
                events_this_year = events.get_event_data(
                    ['event_number', 'event_year', 'event_magnitude'],
                    event_range = {'type':'year', 'filter':(current_year - 1, current_year)}
                )

            progress_indicator_year = current_year - 1 - start_year + year_frame/fpy
            
            evnums_this_frame = []
            sids_this_frame = []
            for i, year in enumerate(events_this_year['event_year']):
                if math.modf(year)[0] <= float(year_frame+1)/fpy and math.modf(year)[0] > float(year_frame)/fpy:
                    evnums_this_frame.append(events_this_year['event_number'][i])
                    sids_this_frame.append(geometry.sections_with_elements(list(events.get_event_elements(events_this_year['event_number'][i]))))
                    cumulative_magnitudes.append(events_this_year['event_magnitude'][i])
            sids_this_frame = set( itertools.chain(*sids_this_frame) )

            sys.stdout.write('frame {} (year {}) of {} ({})\n'.format(the_frame, progress_indicator_year, total_frames, total_years))

            # Remove a fixed percentage from the field. This is the decay
            # that slowly fades existing field values.
            EF.shrink_field(1.0 - 1.0/(animation_fps*fade_seconds))
            
            #-------------------------------------------------------------------
            # Load or calculate all of the data for the current frame.
            #-------------------------------------------------------------------
            if len(evnums_this_frame) > 0:
                for i, evnum in enumerate(evnums_this_frame):
                    sys.stdout.write('\r Event {} :: '.format(evnum))
                    # Try and load the fields
                    field_values_loaded = EF.load_field_values('{}{}_'.format(field_values_directory, evnum))
                    if field_values_loaded:
                        sys.stdout.write('loaded'.format(evnum))
                    # If they havent been saved then we need to calculate them
                    elif not field_values_loaded:
                        sys.stdout.write('processing '.format(evnum))
                        sys.stdout.flush()
                        
                        event_element_slips = events.get_event_element_slips(evnum)
                        ele_getter = itemgetter(*event_element_slips.keys())
                        event_element_data = ele_getter(geometry)
                        if len(event_element_slips) == 1:
                            event_element_data = [event_element_data]
                    
                        sys.stdout.write('{} elements :: '.format(len(event_element_slips)))
                        sys.stdout.flush()
                        
                        EF.calculate_field_values(
                            event_element_data,
                            event_element_slips,
                            cutoff=cutoff,
                            save_file_prefix='{}{}_'.format(field_values_directory, evnum)
                        )
                        
                        if i < len(evnums_this_frame)-1 :
                            sys.stdout.write('\033[2K')
                        sys.stdout.flush()

            sys.stdout.write('\n')
            
            #-------------------------------------------------------------------
            # State variables that need to be maintained pre-plotting.
            #-------------------------------------------------------------------
            
            # Big magnitude event label fade state
            for elme in event_large_magnitude_evnums:
                if elme in evnums_this_frame:
                    large_event_label_states[elme] = 1.0
            
            # Active section line width
            for sid in sids_this_frame:
                section_states[sid] = 1.0
            
            #-------------------------------------------------------------------
            # If the image has not been plotted, plot it.
            #-------------------------------------------------------------------
            if not os.path.isfile('{}{}.png'.format(frame_images_directory, the_frame)) or force_plot:
                sys.stdout.write('\r Plotting :: ')
                
                # Set the plot field to be the current field
                EFP.set_field(EF)
                
                sys.stdout.write('creating the map image : ')
                sys.stdout.flush()
                # Create the field map image
                map_image = EFP.create_field_image(fringes=fringes)
                
                sys.stdout.write('finishing plot')
                sys.stdout.flush()
                #---------------------------------------------------------------
                # Plot all of the geographic data, animation timeline etc.
                #---------------------------------------------------------------
                
                # Width and height are fixed
                ph = 768.0
                pw = 1024.0
                
                # Map margin left and right
                mml = 70.0
                mmb = 70.0
                
                # Progress indicator map margin
                pimm = 70.0
                # Progress indicator margin r
                pimr = 50.0
                
                mw = EF.lons_1d.size
                mh = EF.lats_1d.size
                
                width_frac = mw/pw
                height_frac = mh/ph
                left_frac = mml/pw
                bottom_frac = mmb/ph

                pwi = pw/plot_resolution
                phi = ph/plot_resolution
                fig4 = mplt.figure(figsize=(pwi, phi), dpi=plot_resolution)

                #---------------------------------------------------------------
                # m4, fig4 is all of the boundary data.
                #---------------------------------------------------------------
                m4 = Basemap(
                    llcrnrlon=EF.min_lon,
                    llcrnrlat=EF.min_lat,
                    urcrnrlon=EF.max_lon,
                    urcrnrlat=EF.max_lat,
                    lat_0=(EF.max_lat+EF.min_lat)/2.0,
                    lon_0=(EF.max_lon+EF.min_lon)/2.0,
                    resolution=map_resolution,
                    projection=map_projection,
                    suppress_ticks=True
                )
                m4.ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                
                # Draw a frame around the map.
                m4.drawmapboundary(color=map_frame_color, linewidth=map_frame_width, fill_color=(1,1,1,0))

                # draw coastlines, edge of map.
                m4.drawcoastlines(color=coastline_color, linewidth=coastline_width)

                # draw countries
                m4.drawcountries(linewidth=country_width, color=country_color)

                # draw states
                m4.drawstates(linewidth=state_width, color=state_color)

                # draw parallels.
                parallels = np.linspace(EFP.lats_1d.min(), EFP.lats_1d.max(), num_grid_lines+1)
                m4_parallels = m4.drawparallels(
                    parallels,
                    labels=[1,0,0,0],
                    fontsize=map_fontsize,
                    color=grid_color,
                    fontproperties=font,
                    fmt='%.2f',
                    linewidth=grid_width,
                    dashes=[1, 10]
                )

                # draw meridians
                meridians = np.linspace(EFP.lons_1d.min(), EFP.lons_1d.max(), num_grid_lines+1)
                m4_meridians = m4.drawmeridians(
                    meridians,
                    labels=[0,0,1,0],
                    fontsize=map_fontsize,
                    color=grid_color,
                    fontproperties=font,
                    fmt='%.2f',
                    linewidth=grid_width,
                    dashes=[1, 10]
                )

                # add the displacement map image to the plot
                m4.imshow(map_image, origin='upper')
                
                #---------------------------------------------------------------
                # Plot the magnitude/progress indicator.
                #---------------------------------------------------------------
                width_frac = (pw - mw - mml - pimm - pimr) /pw
                height_frac = mh/ph
                left_frac = (mw + mml + pimm)/pw
                bottom_frac = mmb/ph
                
                mag_vs_year = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                mag_vs_year.plot(event_magnitudes, _event_years, color=progress_line_color, linewidth=progress_line_width, marker=progress_marker, mew=progress_marker_edge_width, ms=progress_marker_size, mfc=progress_line_color)
                mag_vs_year.autoscale(enable=True, axis='both', tight=True)
                mag_vs_year.set_ylim(math.floor(min(_event_years)), math.ceil(max(_event_years)))
                mag_vs_year_alt = mag_vs_year.twinx()
                mag_vs_year.set_xticks(np.linspace(min(event_magnitudes),max(event_magnitudes),3))
                
                mag_vs_year_alt.set_ylim(math.floor(min(_event_years)), math.ceil(max(_event_years)))
                
                mag_vs_year_tick_labels = ['{:0.1f}'.format(tick) for tick in np.linspace(min(event_magnitudes),max(event_magnitudes),3)]
                mag_vs_year.set_xticklabels(mag_vs_year_tick_labels)
                
                for tick in mag_vs_year.xaxis.get_major_ticks():
                    tick.label1.set_fontproperties(font)
                    tick.label1.set_fontsize(progress_indicator_fontsize)
                    tick.label1.set_color(progress_tick_color)
                
                for tl in mag_vs_year_alt.get_yticklabels():
                    tl.set_fontproperties(font)
                    tl.set_fontsize(progress_indicator_fontsize)
                    tl.set_color(progress_tick_color)
                
                for tick in mag_vs_year.yaxis.get_major_ticks():
                    tick.label1.set_alpha(0)
                
                for line in mag_vs_year.xaxis.get_ticklines() + mag_vs_year_alt.yaxis.get_ticklines() + mag_vs_year.yaxis.get_ticklines():
                    line.set_alpha(0)

                #take care of all of the frame line widths
                for spine in mag_vs_year.spines.itervalues():
                    spine.set_lw(progress_frame_width)
                    spine.set_color(progress_frame_color)
                
                mag_vs_year.set_xlabel('magnitude', fontproperties=font, size=progress_indicator_fontsize, color=progress_tick_color)
                mag_vs_year_alt.set_ylabel('year', fontproperties=font, size=progress_indicator_fontsize, color=progress_tick_color)

                #---------------------------------------------------------------
                # Add the progress indicator line
                #---------------------------------------------------------------
                mag_vs_year.axhline(y=progress_indicator_year, lw=progress_indicator_linewidth, c=progress_indicator_line_color)
            
                #---------------------------------------------------------------
                # Add the progress indicator label lines for large events.
                #---------------------------------------------------------------
                label_lines = []
                #if len(event_large_magnitudes) < 10:
                #    label_range = _event_large_magnitude_years
                #else:
                label_range = np.linspace(math.floor(min(_event_years)),math.ceil(max(_event_years)),len(event_large_magnitudes))
                for i, y in enumerate(label_range):
                    m = event_large_magnitudes[i]
                    y1 = _event_large_magnitude_years[i]
                    
                    try:
                        lels = large_event_label_states[event_large_magnitude_evnums[i]]
                    except KeyError:
                        lels = 0

                    the_color = mag_color_map(lels)
                    the_line_color = mag_line_colormap(lels)
                    the_line_width = mag_linewidth + lels * (progress_indicator_linewidth)
                    the_bbox = dict(
                        facecolor=current_mag_face_colormap(lels),
                        linewidth=the_line_width,
                        ec=the_line_color,
                        boxstyle='round4,pad=0.5'
                    )
                    mag_vs_year.text(
                        min(event_magnitudes) - 0.4, y, '{:0.1f}'.format(m),
                        fontproperties=font,
                        fontsize=mag_fontsize,
                        horizontalalignment='right',
                        verticalalignment='center',
                        color=the_color,
                        bbox=the_bbox
                    )
                    
                    label_lines.append(
                        mlines.Line2D(
                            [1.0,-0.01,-0.05,-0.085],
                            [y1/math.ceil(max(_event_years)),
                                y1/math.ceil(max(_event_years)),
                                y/math.ceil(max(_event_years)),
                                y/math.ceil(max(_event_years))
                            ],
                            linewidth=the_line_width,
                            transform=mag_vs_year.transAxes,
                            color=the_line_color,
                            solid_capstyle='round',
                            solid_joinstyle='round'
                        )
                    )
                
                mag_vs_year.lines.extend(label_lines)
                
                #---------------------------------------------------------------
                # Plot the current frequency-magnitude
                #---------------------------------------------------------------
                width_frac = 150.0/pw
                height_frac = 150.0/ph
                left_frac = 100.0/pw
                bottom_frac = 100.0/ph
                
                current_total_events = len(cumulative_magnitudes)
                if current_total_events > 1:
                    cum_freq = {}
                    fm_x = []
                    fm_y = []
                    for num, magnitude in enumerate(sorted(cumulative_magnitudes)):
                        cum_freq['{:0.10f}'.format(magnitude)] = current_total_events - (num + 1)
                    
                    for magnitude in sorted(cum_freq.iterkeys()):
                        fm_x.append(magnitude)
                        fm_y.append(float(cum_freq[magnitude]))
                    mag_vs_freq = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                    mag_vs_freq.semilogy(fm_x, fm_y, color=fm_line_color, linewidth=fm_line_width)
                    mag_vs_freq.set_ylim(bottom=min(fm_y_ticks), top=max(fm_y_ticks))
                    mag_vs_freq.set_xlim(left=round(min(event_magnitudes)), right=round(max(event_magnitudes)))
                    mag_vs_freq.set_yticks(fm_y_ticks)
                    mag_vs_freq.set_xticks(fm_x_ticks)
                    
                    mag_vs_freq_x_ticklabels = ['{:0.1f}'.format(tick) for tick in fm_x_ticks]
                    mag_vs_freq.set_xticklabels(mag_vs_freq_x_ticklabels)
                    
                    for tick in mag_vs_freq.xaxis.get_major_ticks() + mag_vs_freq.yaxis.get_major_ticks():
                        tick.label1.set_fontproperties(font)
                        tick.label1.set_fontsize(fm_fontsize)
                        tick.label1.set_color(fm_label_color)
                        tick.label1.set_alpha(fm_alpha_state)
                    
                    for line in mag_vs_freq.xaxis.get_majorticklines() + mag_vs_freq.yaxis.get_majorticklines() + mag_vs_freq.yaxis.get_minorticklines():
                        line.set_alpha(0)
                    
                    #take care of all of the frame line widths
                    for spine in mag_vs_freq.spines.itervalues():
                        spine.set_lw(fm_frame_width)
                        spine.set_color(fm_frame_color)
                        spine.set_alpha(fm_alpha_state)
                    mag_vs_freq.patch.set_alpha(fm_alpha_state)

                #---------------------------------------------------------------
                # Plot the fault traces.
                #---------------------------------------------------------------
                # Plot the sections
                for sid, sec_trace in fault_traces_latlon.iteritems():
                    trace_Xs, trace_Ys = m4(sec_trace[1], sec_trace[0])
                    
                    try:
                        section_state = section_states[sid]
                    except KeyError:
                        section_state = 0.0

                    m4.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=fault_width + section_state * 2.0, solid_capstyle='round', solid_joinstyle='round')

                #---------------------------------------------------------------
                # Plot the colorbar
                #---------------------------------------------------------------
                left_frac = 70.0/pw
                bottom_frac = (70.0 - cb_height - cb_margin_t)/ph
                width_frac = mw/pw
                height_frac = cb_height/ph
                
                cb_ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                norm = EFP.norm
                cb = mcolorbar.ColorbarBase(cb_ax, cmap=cmap,
                       norm=norm,
                       orientation='horizontal')
                if field_type == 'displacement':
                    if fringes:
                        cb_title = 'Displacement [m]'
                    else:
                        cb_title = 'Total displacement [m]'

                elif field_type == 'gravity' or field_type == 'dilat_gravity':
                    cb_title = r'Gravity changes [$\mu gal$]'
                    
                # If contour plot, match colorbar ticks to contour values
                if contours is not None:
                    cb.set_ticks(contours)

                cb_ax.set_title(cb_title, fontproperties=font, color=cb_fontcolor, size=cb_fontsize, va='top', ha='left', position=(0,-1.5) )

                for label in cb_ax.xaxis.get_ticklabels():
                    label.set_fontproperties(font)
                    label.set_fontsize(cb_fontsize)
                    label.set_color(cb_fontcolor)
                                      
                for line in cb_ax.xaxis.get_ticklines():
                    line.set_alpha(0)
            
                #---------------------------------------------------------------
                # If the field is a displacement field, draw the look arrows.
                #---------------------------------------------------------------
                if field_type == 'displacement':
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

                    az_arrow_start_x    = 0.5 - (0.8/2.0)*math.sin(EFP.look_azimuth)
                    az_arrow_start_y    = 0.5 - (0.8/2.0)*math.cos(EFP.look_azimuth)
                    az_arrow_dx      = 0.8*math.sin(EFP.look_azimuth)
                    az_arrow_dy      = 0.8*math.cos(EFP.look_azimuth)

                    az_ax.arrow( az_arrow_start_x , az_arrow_start_y, az_arrow_dx, az_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='right', length_includes_head=True, lw=1.0, fc='k' )
                    az_ax.add_line(mlines.Line2D((0.5,0.5), (0.5,0.8), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
                    az_ax.add_patch(mpatches.Arc((0.5,0.5), 0.3, 0.3, theta1=90.0 - EF.convert.rad2deg(EFP.look_azimuth), theta2=90.0, fc='none', lw=1.0, ls='dotted', ec='k'))
                    az_ax.text(1.0, 1.0, 'az = {:0.1f}{}'.format(EF.convert.rad2deg(EFP.look_azimuth),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')

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

                    al_arrow_start_x    = 0.1 + 0.8*math.cos(EFP.look_elevation)
                    al_arrow_start_y    = 0.1 + 0.8*math.sin(EFP.look_elevation)
                    al_arrow_dx      = -0.8*math.cos(EFP.look_elevation)
                    al_arrow_dy      = -0.8*math.sin(EFP.look_elevation)

                    al_ax.arrow( al_arrow_start_x , al_arrow_start_y, al_arrow_dx, al_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='left', length_includes_head=True, lw=1.0, fc='k' )
                    al_ax.add_line(mlines.Line2D((0.1,0.9), (0.1,0.1), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
                    al_ax.add_patch(mpatches.Arc((0.1,0.1), 0.5, 0.5, theta1=0.0, theta2=EF.convert.rad2deg(EFP.look_elevation), fc='none', lw=1.0, ls='dotted', ec='k'))
                    al_ax.text(1.0, 1.0, 'al = {:0.1f}{}'.format(EF.convert.rad2deg(EFP.look_elevation),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')
                #---------------------------------------------------------------
                # If the field is a gravity field change outermost tick labels
                # on colorbar.
                #---------------------------------------------------------------
                else:
                    # Want to change outermost tick labels on colorbar
                    #   from 'VALUE','-VALUE' to '>VALUE' and '<-VALUE'    
                    cb_tick_labs = [item.get_text() for item in cb_ax.get_xticklabels()]
                    cb_tick_labs[0] = '<'+cb_tick_labs[0]
                    cb_tick_labs[-1] = '>'+cb_tick_labs[-1]
                    cb_ax.set_xticklabels(cb_tick_labs)
                
                #---------------------------------------------------------------
                # Save the figure and clear out all matplotlib figures.
                #---------------------------------------------------------------
                fig4.savefig('{}{}.png'.format(frame_images_directory, the_frame), format='png', dpi=plot_resolution)
                
                fig4.clf()
                mplt.close('all')
                gc.collect()

                sys.stdout.write('\033[2K')
                sys.stdout.flush()
            
            #-------------------------------------------------------------------
            # State variables that need to be maintained post-plotting.
            #-------------------------------------------------------------------
            
            # Big magnitude event label fade state
            for k, lels in large_event_label_states.iteritems():
                if lels >= 0.04:
                    large_event_label_states[k] -= 0.04
                else:
                    large_event_label_states[k] = 0.0
            
            # Active section line width
            for sid, section_state in section_states.iteritems():
                if section_state >= 0.04:
                    section_states[sid] -= 0.04
                else:
                    section_states[sid] = 0.0

            # Frequency magnitude plot alpha
            fm_alpha_state += 0.04
            if fm_alpha_state > 1.0:
                fm_alpha_state = 1.0

            sys.stdout.write('\033[1A')
            sys.stdout.write('\033[2K')
            sys.stdout.write('\033[1A\r')
            sys.stdout.write('\033[2K')
        
        #-----------------------------------------------------------------------
        # Create the movie using ffmpeg.
        #-----------------------------------------------------------------------
        
        #ffmpeg -y -r 15 -i f_%d.png -f mp4 -vcodec h264 -pix_fmt yuv420p animation.mp4
        
        proc_args = 'ffmpeg -y -r {fps} -start_number 0 -i {dir}{inc}.png -f mp4 -vcodec h264 -pix_fmt yuv420p {out}animation.mp4'.format(
            fps=int(animation_fps),
            dir=frame_images_directory,
            inc='%d',
            out=output_directory
        )
        proc = subprocess.Popen(proc_args, shell=True)
        proc.wait()

#-------------------------------------------------------------------------------
# plots event fields
#-------------------------------------------------------------------------------
def plot_event_field(sim_file, evnum, output_directory, field_type='displacement', fringes=True, padding=0.08, cutoff=None, tag=None, hi_res=False,lat_lon=False, angles=None, contours=None, lat_range=None, lon_range=None, wavelength=0.03):
    
    sys.stdout.write('Initializing plot :: ')
    sys.stdout.flush()
    
    field_values_directory = '{}field_values/'.format(output_directory)
    
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        
    if not output_directory.endswith('/'):
        output_directory += '/'
        
    if not os.path.exists(field_values_directory):
        os.makedirs(field_values_directory)
        
    PRE = '{}{}_'.format(field_values_directory, evnum)
      
    if field_type is 'dilat_gravity':
        output_file = output_directory+'{}_dg_dilat'.format(evnum)              
    elif field_type is 'gravity':        
        output_file = output_directory+'{}_dg_total'.format(evnum)
    elif field_type is 'displacement':
        output_file = output_directory+'{}_displace'.format(evnum)

        
    if tag is not None:
        output_file += '_'+tag
        
    output_file += '.png'
        
    
    start_time = time.time()
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events = VCEvents(sim_data)
        geometry = VCGeometry(sim_data)
        
        fault_traces = geometry.get_fault_traces()
    
        event_data = events[evnum]
        event_element_slips = events.get_event_element_slips(evnum)
        ele_getter = itemgetter(*event_element_slips.keys())
        event_element_data = ele_getter(geometry)
        if len(event_element_slips) == 1:
            event_element_data = [event_element_data]

        event_sections = geometry.sections_with_elements(event_element_slips.keys())
        base_lat = geometry.base_lat
        base_lon = geometry.base_lon

        if lat_range is None or lon_range is None:
            min_lat = geometry.min_lat
            max_lat = geometry.max_lat
            min_lon = geometry.min_lon
            max_lon = geometry.max_lon
        else:
            min_lat,max_lat = lat_range
            min_lon,max_lon = lon_range
            # Force map bounds from event elements only
            #converter = quakelib.Conversion(base_lat, base_lon)
            #min_x,max_x,min_y,max_y = geometry.bounds_xyz_from_elements(event_element_data)
            #min_lat = converter.convert2LatLon(quakelib.Vec3(min_x,min_y,0)).lat()
            #max_lat = converter.convert2LatLon(quakelib.Vec3(max_x,max_y,0)).lat()
            #min_lon = converter.convert2LatLon(quakelib.Vec3(min_x,min_y,0)).lon()
            #max_lon = converter.convert2LatLon(quakelib.Vec3(max_x,max_y,0)).lon()

    
        sys.stdout.write('{} elements in {} sections : '.format(len(event_element_slips), len(event_sections)))
        sys.stdout.flush()
        
        sys.stdout.write( '{} field : '.format(field_type))
        sys.stdout.flush()
        
        if field_type == 'displacement':
            EF = vcutils.VCDisplacementField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
        elif field_type == 'gravity':
            EF = vcutils.VCGravityField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding,free_air=True)
        elif field_type == 'dilat_gravity':
            EF = vcutils.VCGravityField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding,free_air=False)


        sys.stdout.write('done\n')
        sys.stdout.flush()

        field_values_loaded = EF.load_field_values(PRE)
        
        if lat_lon:
            EF.save_lat_lon_values(PRE)
            sys.stdout.write('\n***Latitude and longitude values written to file.')
        
     
        if field_values_loaded:
            sys.stdout.write('Loading event {} {} field :: '.format(evnum, field_type))
            sys.stdout.flush()
        if not field_values_loaded:
            sys.stdout.write('Processing event {} {} field :: '.format(evnum, field_type))
            sys.stdout.flush()          
            sys.stdout.write('{} elements : '.format(len(event_element_slips)))
            sys.stdout.flush()
            
            EF.calculate_field_values(
                        event_element_data,
                        event_element_slips,
                        cutoff=cutoff,
                        save_file_prefix=PRE)

        
        sys.stdout.write('done\n')
        sys.stdout.flush()

        sys.stdout.write('Plotting :: initializing : ')
        sys.stdout.flush()

        if field_type == 'displacement':
            EFP = vcplotutils.VCDisplacementFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon, wavelength=wavelength)
            if angles is None:
                EFP.calculate_look_angles(geometry[:])
            else:
                EFP.set_look_angles(angles[0], angles[1])
        elif field_type == 'gravity' or field_type == 'dilat_gravity':
            EFP = vcplotutils.VCGravityFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon, contours=contours)

    generate_map(EF,EFP,fault_traces,fringes,event_data,output_file,event_sections,field_type=field_type,hi_res=hi_res)

    """
    EFP.set_field(EF)

    if field_type == 'displacement':
        start_time = time.time()
        EFP.calculate_look_angles(event_element_data)

    sys.stdout.write('map image : ')
    sys.stdout.flush()

    map_image = EFP.create_field_image(fringes=fringes)

    sys.stdout.write('map overlay : ')
    sys.stdout.flush()
    # Convert the fault traces to lat-lon
    fault_traces_latlon = {}
    for secid in fault_traces.iterkeys():
         fault_traces_latlon[secid] = zip(*[(lambda y: (y.lat(),y.lon()))(EF.convert.convert2LatLon(quakelib.Vec3(x[0], x[1], x[2]))) for x in fault_traces[secid]])

    #---------------------------------------------------------------------------
    # Plot all of the geographic info on top of the displacement map image.
    #---------------------------------------------------------------------------
    
    # Grab all of the plot properties that we will need.
    # properties that are fringes dependent
    if fringes and field_type == 'displacement':
        cmap            = EFP.dmc['cmap_f']
        coastline_color = EFP.dmc['coastline_color_f']
        country_color   = EFP.dmc['country_color_f']
        state_color     = EFP.dmc['state_color_f']
        fault_color     = EFP.dmc['fault_color_f']
        map_tick_color  = EFP.dmc['map_tick_color_f']
        map_frame_color = EFP.dmc['map_frame_color_f']
        grid_color      = EFP.dmc['grid_color_f']
        cb_fontcolor    = EFP.dmc['cb_fontcolor_f']
    else:
        cmap            = EFP.dmc['cmap']
        coastline_color = EFP.dmc['coastline_color']
        country_color   = EFP.dmc['country_color']
        state_color     = EFP.dmc['state_color']
        fault_color     = EFP.dmc['fault_color']
        map_tick_color  = EFP.dmc['map_tick_color']
        map_frame_color = EFP.dmc['map_frame_color']
        grid_color      = EFP.dmc['grid_color']
        cb_fontcolor    = EFP.dmc['cb_fontcolor']
    
    # properties that are not fringes dependent
    boundary_width  = EFP.dmc['boundary_width']
    coastline_width = EFP.dmc['coastline_width']
    country_width   = EFP.dmc['country_width']
    state_width     = EFP.dmc['state_width']
    river_width     = EFP.dmc['river_width']
    fault_width     = EFP.dmc['fault_width']
    map_frame_width = EFP.dmc['map_frame_width']
    map_fontsize    = EFP.dmc['map_fontsize']
    arrow_inset     = EFP.dmc['arrow_inset']
    arrow_fontsize  = EFP.dmc['arrow_fontsize']
    cb_fontsize     = EFP.dmc['cb_fontsize']
    cb_height       = EFP.dmc['cb_height']
    cb_margin_t     = EFP.dmc['cb_margin_t']
    grid_width      = EFP.dmc['grid_width']
    num_grid_lines  = EFP.dmc['num_grid_lines']
    font            = EFP.dmc['font']
    font_bold       = EFP.dmc['font_bold']

    map_resolution  = EFP.dmc['map_resolution']
    map_projection  = EFP.dmc['map_projection']
    plot_resolution = EFP.dmc['plot_resolution']

    # The sizing for the image is tricky. The aspect ratio of the plot is fixed,
    # so we cant set all of margins to whatever we want. We will set the anchor
    # to the top, left margin position. Then scale the image based on the
    # bottom/right margin, whichever is bigger.
    
    mw = EF.lons_1d.size
    mh = EF.lats_1d.size

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
        llcrnrlon=EF.min_lon,
        llcrnrlat=EF.min_lat,
        urcrnrlon=EF.max_lon,
        urcrnrlat=EF.max_lat,
        lat_0=(EF.max_lat+EF.min_lat)/2.0,
        lon_0=(EF.max_lon+EF.min_lon)/2.0,
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
    parallels = np.linspace(EFP.lats_1d.min(), EFP.lats_1d.max(), num_grid_lines+1)
    m4_parallels = m4.drawparallels(parallels, labels=[1,0,0,0], fontsize=map_fontsize, color=grid_color, fontproperties=font, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])
    
    # draw meridians
    meridians = np.linspace(EFP.lons_1d.min(), EFP.lons_1d.max(), num_grid_lines+1)
    m4_meridians = m4.drawmeridians(meridians, labels=[0,0,1,0], fontsize=map_fontsize, color=grid_color, fontproperties=font, fmt='%.2f', linewidth=grid_width, dashes=[1, 10])

    if field_type == 'displacement':
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

        az_arrow_start_x    = 0.5 - (0.8/2.0)*math.sin(EFP.look_azimuth)
        az_arrow_start_y    = 0.5 - (0.8/2.0)*math.cos(EFP.look_azimuth)
        az_arrow_dx      = 0.8*math.sin(EFP.look_azimuth)
        az_arrow_dy      = 0.8*math.cos(EFP.look_azimuth)

        az_ax.arrow( az_arrow_start_x , az_arrow_start_y, az_arrow_dx, az_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='right', length_includes_head=True, lw=1.0, fc='k' )
        az_ax.add_line(mlines.Line2D((0.5,0.5), (0.5,0.8), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
        az_ax.add_patch(mpatches.Arc((0.5,0.5), 0.3, 0.3, theta1=90.0 - EF.convert.rad2deg(EFP.look_azimuth), theta2=90.0, fc='none', lw=1.0, ls='dotted', ec='k'))
        az_ax.text(1.0, 1.0, 'az = {:0.1f}{}'.format(EF.convert.rad2deg(EFP.look_azimuth),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')

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

        al_arrow_start_x    = 0.1 + 0.8*math.cos(EFP.look_elevation)
        al_arrow_start_y    = 0.1 + 0.8*math.sin(EFP.look_elevation)
        al_arrow_dx      = -0.8*math.cos(EFP.look_elevation)
        al_arrow_dy      = -0.8*math.sin(EFP.look_elevation)

        al_ax.arrow( al_arrow_start_x , al_arrow_start_y, al_arrow_dx, al_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='left', length_includes_head=True, lw=1.0, fc='k' )
        al_ax.add_line(mlines.Line2D((0.1,0.9), (0.1,0.1), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
        al_ax.add_patch(mpatches.Arc((0.1,0.1), 0.5, 0.5, theta1=0.0, theta2=EF.convert.rad2deg(EFP.look_elevation), fc='none', lw=1.0, ls='dotted', ec='k'))
        al_ax.text(1.0, 1.0, 'al = {:0.1f}{}'.format(EF.convert.rad2deg(EFP.look_elevation),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')
        
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
    for sid, sec_trace in fault_traces_latlon.iteritems():
        trace_Xs, trace_Ys = m4(sec_trace[1], sec_trace[0])
        
        if sid in event_sections:
            linewidth = fault_width + 3
        else:
            linewidth = fault_width

        m4.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=linewidth, solid_capstyle='round', solid_joinstyle='round')

    #plot the cb
    left_frac = 70.0/pw
    bottom_frac = (70.0 - cb_height - cb_margin_t)/ph
    width_frac = mw/pw
    height_frac = cb_height/ph
    
    cb_ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
    norm = EFP.norm
    cb = mcolorbar.ColorbarBase(cb_ax, cmap=cmap,
           norm=norm,
           orientation='horizontal')
    if field_type == 'displacement':
        if fringes:
            cb_title = 'Displacement [m]'
        else:
            cb_title = 'Total displacement [m]'

    elif field_type == 'gravity':
        cb_title        = r'Gravity changes [$\mu gal$]'
        # Make first and last ticks on colorbar be <MIN and >MAX
        cb_tick_labs    = [item.get_text() for item in cb_ax.get_xticklabels()]
        cb_tick_labs[0] = '<'+cb_tick_labs[0]
        cb_tick_labs[-1]= '>'+cb_tick_labs[-1]
        cb_ax.set_xticklabels(cb_tick_labs)

    cb_ax.set_title(cb_title, fontproperties=font, color=cb_fontcolor, size=cb_fontsize, va='top', ha='left', position=(0,-1.5) )

    for label in cb_ax.xaxis.get_ticklabels():
        label.set_fontproperties(font)
        label.set_fontsize(cb_fontsize)
        label.set_color(cb_fontcolor)
    for line in cb_ax.xaxis.get_ticklines():
        line.set_alpha(0)

    if output_file is not None:
        # save the figure
        fig4.savefig(output_file, format='png', dpi=plot_resolution)

    sys.stdout.write('done\n')
    sys.stdout.flush()
    """

#-------------------------------------------------------------------------------
# plots recurrence intervals
#-------------------------------------------------------------------------------
def plot_recurrence_intervals(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None):
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
            section_events = events.get_event_data_from_evnums(
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

    if output_file is not None:
        # Get the plot format and save the file
        plot_format = output_file.split('.')[-1]
        if plot_format != 'png' and plot_format != 'pdf':
            raise vcexceptions.PlotFormatNotSupported(plot_format)
        else:
            fig.savefig(output_file, format=plot_format, dpi=res)

#-------------------------------------------------------------------------------
# plots a matrix of a bipartite event sequence graph
#-------------------------------------------------------------------------------
def plot_bipartite_graph_matrix(graph_file, output_file=None):
    G = cPickle.load(open(graph_file, 'rb'))
    
    top_nodes = set(n for n,d in G.nodes(data=True) if d['bipartite']==0)
    bottom_nodes = [n for n in sorted(set(G) - top_nodes, key=lambda n: G.degree(n), reverse=True)]
            
    print 'calculating matrix'
    bij_matrix = bipartite.biadjacency_matrix(G, bottom_nodes, top_nodes, weight='weight')
    
    print bij_matrix.shape
    # Aspect rations are always h/w
    aspect = float(bij_matrix.shape[1])/float(bij_matrix.shape[0])
    
    print 'plotting matrix', len(top_nodes), len(bottom_nodes), bij_matrix.max()
    #print pos_sid_prob
    
    # plot parameters
    imh = 178.0        # the full image height
    imw = round(imh/aspect)    # the width is set based on the aspect ratio
    lm = 10.0
    rm = 10.0
    tm = 10.0
    bm = 20.0
    res = 72.0
    ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=6)
    cmap = mplt.get_cmap('hot_r')
    
    #arial14 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=14)
    #arial12 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=12)
    #arial10 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=10)
    #arial7_light = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=7, weight='light')
    
    imwi = imw/res
    imhi = imh/res
    fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
    ph = imh - tm - bm # the height for both matricies
    pw = imw - lm - rm
    print imh, imw, lm/imw, bm/imh, pw/imw, ph/imh
    ax = fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))
    
    #print pos_sid_prob[::-1]
    
    ax.imshow(bij_matrix.T, interpolation='none', cmap=cmap)
    #ax.set_xticks(range(len(pos_sid_prob)))
    #ax.set_yticks(range(len(pos_sid_prob)))
    #ax.set_xticklabels(pos_sid_prob)
    #ax.set_yticklabels(pos_sid_prob)
    #for label in ax.xaxis.get_ticklabels()+ax.yaxis.get_ticklabels():
    #    label.set_fontproperties(ticklabelfont)

    #ax.axis('tight')
    #ax.set_ylim((15.5, -0.5))
    #ax.set_xlim((140.5, len(pos_sid_prob)-0.5))

    if output_file is not None:
        # save the figure
        fig.savefig(output_file, format='png', dpi=res)

#-------------------------------------------------------------------------------
# plots a matrix of an event graph
#-------------------------------------------------------------------------------
def plot_graph_matrix(graph_file, output_file=None):
    G = cPickle.load(open(graph_file, 'rb'))
    
    matrix_prob, pos_sid_prob = nx.attr_matrix(G, edge_attr='weight', normalized=False)
    #matrix_mean, pos_sid_mean = nx.attr_matrix(G, edge_attr='duration_mean')
    #matrix_std, pos_sid_std = nx.attr_matrix(G, edge_attr='duration_std')
    '''
    for n in nx.nodes_iter(G):
        if G.node[n]['type'] == 'section':
            sequences_by_degree[n] = G.degree(n)
            
    print 'calculating matrix'
    bij_matrix = nx.attr_matrix(G, rc_order=['section', 'sequence'])
    
    print np.max(matrix_prob)
    '''
    '''
    print 'plotting matrix', matrix_prob.shape
    #print pos_sid_prob
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
    ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=6)
    
    #arial14 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=14)
    #arial12 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=12)
    #arial10 = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=10)
    #arial7_light = mpl.font_manager.FontProperties(family='Arial', style='normal', variant='normal', size=7, weight='light')
    
    imwi = imw/res
    imhi = imh/res
    fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
    ph = imh - tm - bm - cbh - cbs
    pw = imw - lm - rm
    ax = fig.add_axes((lm/imw, (bm+cbh+cbs)/imh, pw/imw, ph/imh))
    
    #print pos_sid_prob[::-1]
    
    ax.imshow(matrix_prob[(0,1,2,3,4,5,6,7,8),:], aspect='auto', interpolation='none', cmap = mplt.get_cmap('hot_r'))
    '''
    ax.set_xticks(range(13))
    ax.set_yticks(range(13))
    ax.set_xticklabels(pos_sid_prob[0:13])
    ax.set_yticklabels(pos_sid_prob[0:13])
    ax.set_ylim((12.5, -0.5))
    ax.set_xlim((-0.5, 12.5))
    '''
    #fig.savefig('local/pcolor_test.pdf', format='pdf', dpi=res)
    
    '''
    for sid in pos_sid_prob:
        if G.has_edge(sid,13):
            print sid, G.edge[sid][13]['weight']
        else:
            print sid, 0
    '''
    #ax.set_xticks(range(len(pos_sid_prob)))
    #ax.set_yticks(range(len(pos_sid_prob)))
    #ax.set_xticklabels(pos_sid_prob)
    #ax.set_yticklabels(pos_sid_prob)
    #for label in ax.xaxis.get_ticklabels()+ax.yaxis.get_ticklabels():
    #    label.set_fontproperties(ticklabelfont)

    #ax.axis('tight')
    #ax.set_ylim((15.5, -0.5))
    #ax.set_xlim((140.5, len(pos_sid_prob)-0.5))

#-------------------------------------------------------------------------------
# plots an event graph
#-------------------------------------------------------------------------------
def plot_graph(graph_file, output_file, degree_cut=None, label_degree_cut=0.25, self_loops=True, sequence=None):
    if type(graph_file) is str:
        G = cPickle.load(open(graph_file, 'rb'))
    else:
        G = graph_file
    
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
    node_min = 0.05
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
        try:
            sequence_node = Gsub.node[n]['bipartite']
        except KeyError:
            sequence_node = False

        if not sequence_node:
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
        else:
            widths[n] = node_min
            heights[n] = node_min
            styles[n] = 'filled'
            colors[n] = '#FFFFF'
            node_line_widths[n] = 0.5
            labels[n] = ''

    nx.set_node_attributes(Gsub,'width',widths)
    nx.set_node_attributes(Gsub,'height',heights)
    nx.set_node_attributes(Gsub,'label',labels)
    nx.set_node_attributes(Gsub,'style',styles)
    nx.set_node_attributes(Gsub,'fillcolor',colors)
    nx.set_node_attributes(Gsub,'penwidth',node_line_widths)
    nx.set_node_attributes(Gsub,'fontsize',font_sizes)
    #print G.edges(data=True)
    
    try:
        weights = [ float(edata['weight']) for u,v,edata in Gsub.edges(data=True) ]
    except KeyError:
        weights = [ 1.0 for u,v in Gsub.edges() ]
    
    max_weight = float(max(weights))
    min_weight = float(min(weights))
    line_min = 0.1
    line_max = 5.0
    arrow_min = 0.1
    arrow_max = 0.7
    alpha_min = 10.0
    alpha_max = 255.0

    # all of the scales will be at their max
    if min_weight == max_weight:
        min_weight = 0.0
        line_max = 0.5
        arrow_max = 0.1
        alpha_max = 50.0
    
    edge_widths = {}
    arrow_sizes = {}
    edge_colors = {}

    if sequence is not None:
        sequence_parsed = []
        for n,sid in enumerate(sequence):
            if n < len(sequence) - 1:
                sequence_parsed.append('{}-{}'.format(sid, sequence[n+1]))

    for e in nx.edges_iter(Gsub):
        try:
            e_weight = float(Gsub[e[0]][e[1]]['weight'])
        except KeyError:
            e_weight = 1.0

        width = vcutils.linear_interp(e_weight, min_weight, max_weight, line_min, line_max)
        arrow_size = vcutils.linear_interp(e_weight, min_weight, max_weight, arrow_min, arrow_max)
        alpha = vcutils.linear_interp(e_weight, min_weight, max_weight, alpha_min, alpha_max)
        edge_widths[e] = width
        arrow_sizes[e] = arrow_size
        if sequence is not None:
            if '{}-{}'.format(e[0], e[1]) in sequence_parsed:
                edge_colors[e] = '#FF0000'
            else:
                edge_colors[e] = '#000000{:x}'.format(int(alpha))
        else:
            edge_colors[e] = '#000000{:x}'.format(int(alpha))
    
    nx.set_edge_attributes(Gsub, 'penwidth', edge_widths)
    nx.set_edge_attributes(Gsub, 'arrowsize', arrow_sizes)
    nx.set_edge_attributes(Gsub, 'color', edge_colors)
    #cmap = mplt.get_cmap('gray')
    
    #norm = mcolor.Normalize(vmin=min(edge_weights), vmax=max(edge_weights))
    
    A=nx.to_agraph(Gsub)        # convert to a graphviz graph
    #A.layout()            # neato layout
    A.draw(output_file, prog='sfdp', args='-Gsize="40!" -Goverlap=prism -Grepulsiveforce=3.0 -GsmoothType="graph_dist" -Goutputorder="edgesfirst" -Nfixedsize="true" -Nfontname="Helvetica"')

#-------------------------------------------------------------------------------
# space-time plot
#-------------------------------------------------------------------------------
def space_time_plot(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None):
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
        stp = vcplotutils.VCSpaceTimePlot(
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
            '''
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
            '''
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
def magnitude_rupture_area(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None, scale=False):
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
    event_area_kmsq = [quakelib.Conversion().sqm2sqkm(x) for x in event_data['event_area']]
    
    # get the binned averages of the data
    x_ave, y_ave = vcplotutils.calculate_averages(event_area_kmsq, event_data['event_magnitude'])
    
    # get the plot label which will depend on the filters
    plot_label = vcplotutils.get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)
    
    x_WC = np.linspace(2.2,5184)
    y_WC = 4.07 + 0.98 * np.log10(x_WC)
    y_error_plus_WC = 4.07+0.06 + (0.98+0.03) * np.log10(x_WC)
    y_error_minus_WC = 4.07-0.06 + (0.98-0.03) * np.log10(x_WC)
    y_error_WC = [np.subtract(y_WC, y_error_minus_WC), np.subtract(y_error_plus_WC, y_WC)]

    # do the standard plot
    vcplotutils.standard_plot(output_file, event_area_kmsq, event_data['event_magnitude'],
        axis_format='semilogx',
        scale=scale,
        add_lines=[
            {'label':'binned average', 'x':x_ave, 'y':y_ave},
            {'label':'WC', 'x':x_WC, 'y':y_WC, 'ls':'--', 'c':'red'}
        ],
        axis_labels = {'x':r'Rupture Area [km$^\mathsf{2}$]', 'y':'Magnitude'},
        plot_label='Magnitude-Rupture Area{}'.format(plot_label)
    )
        
#-------------------------------------------------------------------------------
# magnitude average slip plot
#-------------------------------------------------------------------------------
def magnitude_average_slip(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None, scale=False):
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
    x_ave, y_ave = vcplotutils.calculate_averages(event_data['event_average_slip'], event_data['event_magnitude'])
    
    # get the plot label which will depend on the filters
    plot_label = vcplotutils.get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)

    x_WC = np.linspace(0.05,8, num=10)
    y_WC = 6.93 + 0.82 * np.log10(x_WC)
    y_error_plus_WC = 6.93+0.05 + (0.82+0.1) * np.log10(x_WC)
    y_error_minus_WC = 6.93-0.05 + (0.82-0.1) * np.log10(x_WC)
    y_error_WC = [np.subtract(y_WC, y_error_minus_WC), np.subtract(y_error_plus_WC, y_WC)]
    
    # do the standard plot
    vcplotutils.standard_plot(output_file, event_data['event_average_slip'], event_data['event_magnitude'],
        axis_format='semilogx',
        scale=scale,
        add_lines=[
            {'label':'binned average', 'x':x_ave, 'y':y_ave},
            {'label':'WC', 'x':x_WC, 'y':y_WC, 'ls':'--', 'c':'red'}
        ],
        axis_labels = {'y':'Magnitude', 'x':'Average Slip [m]'},
        plot_label='Magnitude-Average Slip{}'.format(plot_label)
    )

#-------------------------------------------------------------------------------
# average slip surface rupture length plot
#-------------------------------------------------------------------------------
def average_slip_surface_rupture_length(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None, scale=False):
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
    event_surface_rupture_length_km = [quakelib.Conversion().m2km(x) for x in event_data['event_surface_rupture_length']]
    
    # get the binned averages of the data
    x_ave, y_ave = vcplotutils.calculate_averages(event_surface_rupture_length_km, event_data['event_average_slip'])
    
    # get the plot label which will depend on the filters
    plot_label = vcplotutils.get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)

    x_WC = np.linspace(3.8,432, num=10)
    y_WC = 10.0**(-1.43 + 0.88 * np.log10(x_WC))
    y_error_plus_WC = 10.0**(-1.43+0.18 + (0.88+0.11) * np.log10(x_WC))
    y_error_minus_WC = 10.0**(-1.43-0.18 + (0.88-0.11) * np.log10(x_WC))
    y_error_WC = [np.subtract(y_WC, y_error_minus_WC), np.subtract(y_error_plus_WC, y_WC)]
    
    # do the standard plot
    vcplotutils.standard_plot(output_file, event_surface_rupture_length_km, event_data['event_average_slip'],
        axis_format='loglog',
        scale=scale,
        add_lines=[
            {'label':'binned average', 'x':x_ave, 'y':y_ave},
            {'label':'WC', 'x':x_WC, 'y':y_WC, 'ls':'--', 'c':'red'}
        ],
        axis_labels = {'y':'log(Average Slip [m])', 'x':'Surface Rupture Length [km]'},
        plot_label='Average Slip-Surface Rupture Length{}'.format(plot_label)
    )

#-------------------------------------------------------------------------------
# frequency magnitude plot
#-------------------------------------------------------------------------------
def frequency_magnitude(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None, scale=False):
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
    plot_label = vcplotutils.get_plot_label(sim_file, event_range=event_range, section_filter=section_filter, magnitude_filter=magnitude_filter)
    
    # for the UCERF2 error bars
    x_UCERF = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
    y_UCERF = [4.73, 2.15, 0.71, 0.24, 0.074, 0.020]
    y_error_UCERF = [[1.2, 0.37, 0.22, 0.09, 0.04, 0.016],[1.50, 0.43, 0.28, 0.11, 0.06, 0.035]]
    
    # do the standard plot
    vcplotutils.standard_plot(output_file, x, y,
        axis_format='semilogy',
        scale=scale,
        #add_lines=[{'label':'b=1', 'x':x_b1, 'y':y_b1}, {'label':'UCERF2', 'x':x_UCERF, 'y':y_UCERF, 'ls':'--', 'c':'red', 'y_error':y_error_UCERF}],
        add_lines=[{'label':'UCERF2', 'x':x_UCERF, 'y':y_UCERF, 'ls':'--', 'c':'red', 'y_error':y_error_UCERF}],
        axis_labels = {'y':'events per year', 'x':'Magnitude'},
        plot_label='Frequency-Magnitude{}'.format(plot_label),
        connect_points=True,
        legend_loc='upper right'
    )

#-------------------------------------------------------------------------------
# event return map
#-------------------------------------------------------------------------------
def event_return_map(sim_file, output_file=None, event_range=None, section_filter=None, magnitude_filter=None):
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
        event_data = events.get_event_data(['event_trigger', 'event_year', 'event_magnitude', 'event_number'], event_range=event_range, magnitude_filter=magnitude_filter, section_filter=section_filter)
        
        xs = []
        ys = []
        sections = set()
        mags = []
        max_mag_on_section = {}
        
        
        # add edges and nodes to the graph for each event
        for i, event_trigger in enumerate(event_data['event_trigger']):
            if i%round(float(len(event_data['event_year']))/100.0) == 0:
                sys.stdout.write('\r event {} of {}'.format(i, len(event_data['event_year'])))
                sys.stdout.flush()
            
            this_sid = geometry.sections_with_elements([event_trigger])[0]
            sections.add(this_sid)
            if i < len(event_data['event_trigger']) - 1:
                try:
                    if event_data['event_magnitude'][i] > max_mag_on_section[this_sid]:
                        max_mag_on_section[this_sid] = event_data['event_magnitude'][i]
                except KeyError:
                    max_mag_on_section[this_sid] = event_data['event_magnitude'][i]
                
                mags.append(max_mag_on_section[this_sid])
                
                next_sid = geometry.sections_with_elements([event_data['event_trigger'][i+1]])[0]
                xs.append(this_sid)
                ys.append(next_sid)
    
    labels = []

    #print len(xs), len(ys), len(event_data['event_magnitude'][0:-1])
    
    for index in range(max(sections)):
        if index + 1 in sections:
            labels.append('{}'.format(index+1))
        else:
            labels.append('')
    
    # plot parameters
    imw = 1024.0 # the full image width
    imh = 1024.0
    lm = 40.0
    rm = 50.0
    tm = 50.0
    bm = 50.0
    res = 72.0
    ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=7)
    
    imwi = imw/res
    imhi = imh/res
    fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
    ph = imh - tm - bm # the height for both matricies
    pw = imw - lm - rm
    ax = fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))

    #ax.plot(time, output)
    #ax.set_ylim((1, max(pos_sid)))
    
    norm = mcolor.Normalize(vmin=min(mags), vmax=max(mags))
    cmap = mplt.get_cmap('YlOrRd')
    
    ax.scatter(xs,ys, c=cmap(norm(mags)), marker='s', linewidths=1, edgecolors=cmap(norm(mags)))
    ax.set_xticks(range(max(sections)))
    ax.set_yticks(range(max(sections)))
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels)
    ax.axis('tight')

    for label in ax.xaxis.get_ticklabels()+ax.yaxis.get_ticklabels():
        label.set_fontproperties(ticklabelfont)
    ax.set_ylim((0, max(sections)))
    ax.set_xlim((0, max(sections)))




#-------------------------------------------------------------------------------
# event field evolution
#-------------------------------------------------------------------------------
def event_field_evolution(sim_file, output_directory, sim_time_range,
    field_type='gravity', fringes=True, padding=0.08, cutoff=None,
    animation_target_length=50.0, animation_fps = 20.0,
    min_mag_marker = 6.5, start_year = 0.0, duration=100.0,force_plot=False,section_filter=None):
    
    # ----------------------------- Initializing --------------------------
    sys.stdout.write('Initializing animation :: ')
    sys.stdout.flush()
    
    # create the animation dir if needed
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    if not output_directory.endswith('/'):
        output_directory += '/'
    
    # create the animation subdirs if needed
    field_values_directory = '{}field_values/'.format(output_directory)
    frame_images_directory = '{}frame_images/'.format(output_directory)
    if not os.path.exists(field_values_directory):
        os.makedirs(field_values_directory)
    if not os.path.exists(frame_images_directory):
        os.makedirs(frame_images_directory)

    #---------------------------------------------------------------------------
    # Open the data file. It needs to stay open while we do the animation.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events      = VCEvents(sim_data)
        geometry    = VCGeometry(sim_data)
        
        # Get global information about the simulations geometry
        min_lat     = geometry.min_lat
        max_lat     = geometry.max_lat
        min_lon     = geometry.min_lon
        max_lon     = geometry.max_lon
        base_lat    = geometry.base_lat
        base_lon    = geometry.base_lon
        fault_traces= geometry.get_fault_traces()
        
        # Get event information, filter by section if specified
        event_data = events.get_event_data(['event_magnitude', 'event_year', 'event_number'], event_range=sim_time_range,section_filter=section_filter)


        # Get elements and slips for events
        event_element_slips = {evid:events.get_event_element_slips(evid) for evid in event_data['event_number']}
        
        # The blocks in slip_rates define the elements that are used for the plotting
        slip_rates          = geometry.get_slip_rates(section_filter=section_filter,per_year=True)


        # Extract relevant data from events inthe time range of the animation
        event_magnitudes    = event_data['event_magnitude']
        event_years         = event_data['event_year']
        event_numbers       = event_data['event_number']
        current_year        = start_year 
        dt                  = float(duration)/(float(animation_target_length)*animation_fps)
        
        # These are the event years shifted so the begin at zero.
        _event_years = [y - start_year for y in event_years]
        
        # The large magnitudes and shifted years to be marked on the timeline.
        event_large_magnitudes       = [m for m in event_magnitudes if m > min_mag_marker]
        _event_large_magnitude_years = [_event_years[i_m[0]] for i_m in enumerate(event_magnitudes) if i_m[1] >= min_mag_marker]
        event_large_magnitude_evnums = [event_numbers[i_m[0]] for i_m in enumerate(event_magnitudes) if i_m[1] > min_mag_marker]
        
        # Calculate the frames per year and the total number of frames
        #     total years start at start_year and ends 5 years after the last event
        #N_time_steps = duration/float(dt)
        total_years  = float(duration)
        fpy          = math.ceil(animation_target_length*animation_fps/total_years)
        total_frames = int(fpy * total_years)
        
        
        #if total_frames != N_time_steps:
        #    sys.stdout.write('\nHAD TO SET NUM FRAMES MANUALLY!')
        #    total_frames = N_time_steps
        
        # Instantiate the field and the plotter
        if field_type == 'displacement':
            EF = vcutils.VCDisplacementField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
            EFP = vcplotutils.VCDisplacementFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)
            EFP.calculate_look_angles(geometry[:])
        elif field_type == 'gravity':
            EF = vcutils.VCGravityField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
            EFP = vcplotutils.VCGravityFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)

        # I do not know if the following is necessary for my purposes
        #-----------------------------------------------------------------------
        # Find the biggest event and normalize based on these values.
        #-----------------------------------------------------------------------
        """
        if field_type == 'displacement' and not fringes or field_type == 'gravity':
            sys.stdout.write('normalizing : ')
            sys.stdout.flush()
            max_mag_evnum = event_numbers[event_magnitudes.index(max(event_magnitudes))]
            field_values_loaded = EF.load_field_values('{}{}_'.format(field_values_directory, max_mag_evnum))
            if field_values_loaded:
                sys.stdout.write('event {} loaded : '.format(max_mag_evnum))
                sys.stdout.flush()
            if not field_values_loaded:
                sys.stdout.write('event {} processing : '.format(max_mag_evnum))
                sys.stdout.flush()
                event_element_slips = events.get_event_element_slips(max_mag_evnum)
                ele_getter = itemgetter(*event_element_slips.keys())
                event_element_data = ele_getter(geometry)
                if len(event_element_slips) == 1:
                    event_element_data = [event_element_data]
            
                sys.stdout.write('{} elements : '.format(len(event_element_slips)))
                sys.stdout.flush()
                
                EF.calculate_field_values(
                    event_element_data,
                    event_element_slips,
                    cutoff=cutoff,
                    save_file_prefix='{}{}_'.format(field_values_directory,max_mag_evnum)
                )
            EFP.set_field(EF)
            EFP.create_field_image()
        """
        
        
        # Convert the fault traces to lat-lon
        fault_traces_latlon = {}
        for secid in fault_traces.iterkeys():
             fault_traces_latlon[secid] = zip(*[(lambda y: (y.lat(),y.lon()))(EF.convert.convert2LatLon(quakelib.Vec3(x[0], x[1], x[2]))) for x in fault_traces[secid]])

        # Grab all of the plot properties that we will need.
        # properties that are fringes dependent
        if fringes and field_type == 'displacement':
            cmap            = EFP.dmc['cmap_f']
            coastline_color = EFP.dmc['coastline_color_f']
            country_color   = EFP.dmc['country_color_f']
            state_color     = EFP.dmc['state_color_f']
            fault_color     = EFP.dmc['fault_color_f']
            map_tick_color  = EFP.dmc['map_tick_color_f']
            map_frame_color = EFP.dmc['map_frame_color_f']
            grid_color      = EFP.dmc['grid_color_f']
            cb_fontcolor    = EFP.dmc['cb_fontcolor_f']
        else:
            cmap            = EFP.dmc['cmap']
            coastline_color = EFP.dmc['coastline_color']
            country_color   = EFP.dmc['country_color']
            state_color     = EFP.dmc['state_color']
            fault_color     = EFP.dmc['fault_color']
            map_tick_color  = EFP.dmc['map_tick_color']
            map_frame_color = EFP.dmc['map_frame_color']
            grid_color      = EFP.dmc['grid_color']
            cb_fontcolor    = EFP.dmc['cb_fontcolor']
        
        # properties that are not fringes dependent
        boundary_width  = EFP.dmc['boundary_width']
        coastline_width = EFP.dmc['coastline_width']
        country_width   = EFP.dmc['country_width']
        state_width     = EFP.dmc['state_width']
        river_width     = EFP.dmc['river_width']
        fault_width     = EFP.dmc['fault_width']
        map_frame_width = EFP.dmc['map_frame_width']
        map_fontsize    = EFP.dmc['map_fontsize']
        arrow_inset     = EFP.dmc['arrow_inset']
        arrow_fontsize  = EFP.dmc['arrow_fontsize']
        cb_fontsize     = EFP.dmc['cb_fontsize']
        cb_height       = EFP.dmc['cb_height']
        cb_margin_t     = EFP.dmc['cb_margin_t']
        grid_width      = EFP.dmc['grid_width']
        num_grid_lines  = EFP.dmc['num_grid_lines']
        font            = EFP.dmc['font']
        font_bold       = EFP.dmc['font_bold']

        map_resolution  = EFP.dmc['map_resolution']
        map_projection  = EFP.dmc['map_projection']
        plot_resolution = EFP.dmc['plot_resolution']

        #animation specific properties
        progress_tick_color = 'k'
        progress_frame_color = 'k'
        progress_frame_width = 1
        progress_line_color = 'k'
        progress_line_width = 0
        progress_marker = '.'
        progress_marker_edge_width = 0
        progress_marker_size = 4
        progress_indicator_line_color = 'red'
        progress_indicator_linewidth = 2
        progress_indicator_fontsize = 10.0
        
        mag_color = 'k'
        current_mag_color = 'white'
        current_mag_facecolor = 'red'
        mag_facecolor = 'white'
        mag_linewidth = 0.5
        mag_fontsize = 12
        
        fm_fontsize = 9
        if field_type=='displacement':
            fm_label_color = 'white'
        else:
            fm_label_color = 'black'
            
        fm_frame_width = 1
        fm_frame_color = 'k'
        fm_line_color = '0.0'
        fm_line_width = 1
        
        mag_color_map = mcolor.LinearSegmentedColormap.from_list(
            'mag_color_map',
            [mag_color,current_mag_color],
            N=256,
            gamma=1.0
        )
        mag_line_colormap = mcolor.LinearSegmentedColormap.from_list(
            'mag_line_colormap',
            [progress_frame_color,progress_indicator_line_color],
            N=256,
            gamma=1.0
        )
        current_mag_face_colormap = mcolor.LinearSegmentedColormap.from_list(
            'current_mag_face_colormap',
            [mag_facecolor,current_mag_facecolor],
            N=256,
            gamma=1.0
        )

        # Set up the field values. In order to do smooth fading, we are gonna be
        # adding to these quantities when there are new events and subtracting a
        # bit each frame till we get back down to zero.
        EF.init_field(0.0)

        # We need to keep track of all of the magnitudes that have been plotted
        # so far for the FM plot. Also, set up some other properties of the FM
        # plot.
        cumulative_magnitudes = []
        fm_x_ticks = np.linspace(round(min(event_magnitudes)), round(max(event_magnitudes)), 5)
        fm_y_ticks = np.logspace(0, 3, 4)
        
        # Set up a cuple of dicts to maintain the state of the large magnitude
        # labels and the section lines. This will allow us to fade these items
        # out after an event.
        section_states = {}
        large_event_label_states = {}
        fm_alpha_state = 0.0
        
        sys.stdout.write('done\n')
        sys.stdout.flush()
        
        #-----------------------------------------------------------------------
        # Go through all of the frames.
        #    TREAT EACH FRAME AS THE TIME STEP SEPARATING VALUES IN THE 
        #          SLIP TIME SERIES. EVENT SLIPS ARE INCLUDED IN THE 
        #          SLIP TIME SERIES.
        #-----------------------------------------------------------------------
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        sys.stdout.write('Total frames : {}, Frames per year : {}\n'.format(total_frames, fpy))
        sys.stdout.flush()
        for the_frame in range(total_frames):
            year_frame = the_frame%fpy
            """
            # Evaluate the slips on the right endpoint of each dt interval
            current_year += dt
            
            # Grab element geometry data, for elements specified in slip_rates
            ele_getter         = itemgetter(*slip_rates.keys())
            element_data       = ele_getter(geometry)
            
            # Back-slip the elements. If this is first frame, get the current slip by
            #     summing over seismic history up to present
            if the_frame==0:
                frame_slips       = {bid:-1.0*current_year*float(slip_rates[bid]) for bid in slip_rates.keys()}
                
                events_this_frame = events.get_event_data(
                    ['event_number', 'event_year', 'event_magnitude'],
                    event_range = {'type':'year', 'filter':(0, current_year)})
            else:
                # Back slip by 1 time step
                frame_slips       = {bid:-1.0*dt*float(slip_rates[bid]) for bid in slip_rates.keys()}
            
                # Grab the event element slips for the current frame
                events_this_frame = events.get_event_data(
                    ['event_number', 'event_year', 'event_magnitude'],
                    event_range = {'type':'year', 'filter':(current_year - dt, current_year)})
            
            progress_indicator_year = current_year - start_year
            
            evnums_this_frame = []
            sids_this_frame = []
            
            # Gather information for events in this frame
            if len(events_this_frame.keys()) >= 1:
                for i, year in enumerate(events_this_frame['event_year']):
                    evnums_this_frame.append(events_this_frame['event_number'][i])
                    sids_this_frame.append(geometry.sections_with_elements(list(events.get_event_elements(events_this_frame['event_number'][i]))))
                    cumulative_magnitudes.append(events_this_frame['event_magnitude'][i])

                sids_this_frame = set( itertools.chain(*sids_this_frame) )
            

            sys.stdout.write('frame {} (year {}) of {} ({})\n'.format(the_frame, progress_indicator_year, total_frames, total_years))
            
            #-------------------------------------------------------------------
            # Try and load the fields
            #-------------------------------------------------------------------
            field_values_loaded = EF.load_field_values('{}{}_'.format(field_values_directory, the_frame))
                
            if field_values_loaded:
                sys.stdout.write('loaded '.format(the_frame))
            # If they havent been saved then we need to calculate them
            elif not field_values_loaded:
                sys.stdout.write('processing '.format(the_frame))
                sys.stdout.flush()

                #-------------------------------------------------------------------
                # If there are any events, apply event slips to the elements 
                #-------------------------------------------------------------------
                if len(evnums_this_frame) >= 1:
                    # Loop over the events in the current frame    
                    for evid in evnums_this_frame:
                        # For element in the event, add any co-seismic slips
                        for bid in event_element_slips[evid].keys():
                            frame_slips[bid] += float(event_element_slips[evid][bid])
                        
    
                sys.stdout.write('{} elements :: '.format(len(frame_slips)))
                sys.stdout.flush()

                #---------------------------------------------------------------------
                # Now pass the accumulated slips to the greens functions for plotting
                #---------------------------------------------------------------------
                EF.calculate_field_values(
                                element_data,
                                frame_slips,
                                cutoff=cutoff,
                                save_file_prefix='{}{}_'.format(field_values_directory, the_frame),
                )
            

            sys.stdout.write('\n')
            """
            # The current time, evaluated on left side of each time interval
            frame_time = float(start_year + the_frame*dt)
            
            if year_frame == 0:
                current_year += 1
                events_this_year = events.get_event_data(
                    ['event_number', 'event_year', 'event_magnitude'],
                    event_range = {'type':'year', 'filter':(current_year - 1, current_year)}
                )
                
            #-------------------------------------------------------------------    
            # apply back slip to all elements
            #-------------------------------------------------------------------
            
            frame_slips  = {bid:-1.0*frame_time*float(slip_rates[bid]) for bid in slip_rates.keys()}

            progress_indicator_year = current_year - 1 - start_year + year_frame/fpy
            
            evnums_this_frame = []
            sids_this_frame = []
            for i, year in enumerate(events_this_year['event_year']):
                if math.modf(year)[0] <= float(year_frame+1)/fpy and math.modf(year)[0] > float(year_frame)/fpy:
                    this_evid = events_this_year['event_number'][i]
                    
                    #-------------------------------------------------------
                    # Grab the event element slips and add to the back slip 
                    for bid in event_element_slips[this_evid].keys():
                        frame_slips[bid] += float(event_element_slips[this_evid][bid])
                    
                    
                    evnums_this_frame.append(this_evid)
                    sids_this_frame.append(geometry.sections_with_elements(list(events.get_event_elements(this_evid))))
                    cumulative_magnitudes.append(this_evid)
            sids_this_frame = set( itertools.chain(*sids_this_frame) )

            sys.stdout.write('frame {} (year {}) of {} ({})\n'.format(the_frame, progress_indicator_year, total_frames, total_years))

            
            #-------------------------------------------------------------------
            # Load or calculate all of the data for the current frame.
            #-------------------------------------------------------------------

            # Try and load the fields
            field_values_loaded = EF.load_field_values('{}{}_'.format(field_values_directory, the_frame))
            if field_values_loaded:
                sys.stdout.write('loaded frame {}'.format(the_frame))
            # If they havent been saved then we need to calculate them
            elif not field_values_loaded:
                sys.stdout.write('processing frame {}'.format(str(the_frame)))
                sys.stdout.flush()
                
                ele_getter = itemgetter(*frame_slips.keys())
                frame_element_data = ele_getter(geometry)
            
                sys.stdout.write(', {} elements :: '.format(len(frame_slips.keys())))
                sys.stdout.flush()
                
                EF.calculate_field_values(
                    frame_element_data,
                    frame_slips,
                    cutoff=cutoff,
                    save_file_prefix='{}{}_'.format(field_values_directory, the_frame)
                )
                

            sys.stdout.write('\n')

            #=====================================================================
            
            #-------------------------------------------------------------------
            # State variables that need to be maintained pre-plotting.
            #-------------------------------------------------------------------
            
            # Big magnitude event label 
            for elme in event_large_magnitude_evnums:
                if elme in evnums_this_frame:
                    large_event_label_states[elme] = 1.0
            
            # Active section line width
            for sid in sids_this_frame:
                section_states[sid] = 1.0
            
            #-------------------------------------------------------------------
            # If the image has not been plotted, plot it.
            #-------------------------------------------------------------------
            if not os.path.isfile('{}{}.png'.format(frame_images_directory, the_frame)) or force_plot:
                sys.stdout.write('\r Plotting :: ')
                
                # Set the plot field to be the current field
                EFP.set_field(EF)
                
                sys.stdout.write('creating the map image : ')
                sys.stdout.flush()
                # Create the field map image
                map_image = EFP.create_field_image(fringes=fringes)
                
                sys.stdout.write('finishing plot')
                sys.stdout.flush()
                #---------------------------------------------------------------
                # Plot all of the geographic data, animation timeline etc.
                #---------------------------------------------------------------
                
                # Width and height are fixed
                ph = 768.0
                pw = 1024.0
                
                # Map margin left and right
                mml = 70.0
                mmb = 70.0
                
                # Progress indicator map margin
                pimm = 70.0
                # Progress indicator margin r
                pimr = 50.0
                
                mw = EF.lons_1d.size
                mh = EF.lats_1d.size
                
                width_frac = mw/pw
                height_frac = mh/ph
                left_frac = mml/pw
                bottom_frac = mmb/ph

                pwi = pw/plot_resolution
                phi = ph/plot_resolution
                fig4 = mplt.figure(figsize=(pwi, phi), dpi=plot_resolution)

                #---------------------------------------------------------------
                # m4, fig4 is all of the boundary data.
                #---------------------------------------------------------------
                m4 = Basemap(
                    llcrnrlon=EF.min_lon,
                    llcrnrlat=EF.min_lat,
                    urcrnrlon=EF.max_lon,
                    urcrnrlat=EF.max_lat,
                    lat_0=(EF.max_lat+EF.min_lat)/2.0,
                    lon_0=(EF.max_lon+EF.min_lon)/2.0,
                    resolution=map_resolution,
                    projection=map_projection,
                    suppress_ticks=True
                )
                m4.ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                
                # Draw a frame around the map.
                m4.drawmapboundary(color=map_frame_color, linewidth=map_frame_width, fill_color=(1,1,1,0))

                # draw coastlines, edge of map.
                m4.drawcoastlines(color=coastline_color, linewidth=coastline_width)

                # draw countries
                m4.drawcountries(linewidth=country_width, color=country_color)

                # draw states
                m4.drawstates(linewidth=state_width, color=state_color)

                # draw parallels.
                parallels = np.linspace(EFP.lats_1d.min(), EFP.lats_1d.max(), num_grid_lines+1)
                m4_parallels = m4.drawparallels(
                    parallels,
                    labels=[1,0,0,0],
                    fontsize=map_fontsize,
                    color=grid_color,
                    fontproperties=font,
                    fmt='%.2f',
                    linewidth=grid_width,
                    dashes=[1, 10]
                )

                # draw meridians
                meridians = np.linspace(EFP.lons_1d.min(), EFP.lons_1d.max(), num_grid_lines+1)
                m4_meridians = m4.drawmeridians(
                    meridians,
                    labels=[0,0,1,0],
                    fontsize=map_fontsize,
                    color=grid_color,
                    fontproperties=font,
                    fmt='%.2f',
                    linewidth=grid_width,
                    dashes=[1, 10]
                )

                # add the displacement map image to the plot
                m4.imshow(map_image, origin='upper')
                
                #---------------------------------------------------------------
                # Plot the magnitude/progress indicator.
                #---------------------------------------------------------------
                width_frac = (pw - mw - mml - pimm - pimr) /pw
                height_frac = mh/ph
                left_frac = (mw + mml + pimm)/pw
                bottom_frac = mmb/ph
                
                mag_vs_year = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                mag_vs_year.plot(event_magnitudes, _event_years, color=progress_line_color, linewidth=progress_line_width, marker=progress_marker, mew=progress_marker_edge_width, ms=progress_marker_size, mfc=progress_line_color)
                mag_vs_year.autoscale(enable=True, axis='both', tight=True)
                mag_vs_year.set_ylim(math.floor(min(_event_years)), math.ceil(max(_event_years)))
                mag_vs_year_alt = mag_vs_year.twinx()
                mag_vs_year.set_xticks(np.linspace(min(event_magnitudes),max(event_magnitudes),3))
                
                mag_vs_year_alt.set_ylim(math.floor(min(_event_years)), math.ceil(max(_event_years)))
                
                mag_vs_year_tick_labels = ['{:0.1f}'.format(tick) for tick in np.linspace(min(event_magnitudes),max(event_magnitudes),3)]
                mag_vs_year.set_xticklabels(mag_vs_year_tick_labels)
                
                for tick in mag_vs_year.xaxis.get_major_ticks():
                    tick.label1.set_fontproperties(font)
                    tick.label1.set_fontsize(progress_indicator_fontsize)
                    tick.label1.set_color(progress_tick_color)
                
                for tl in mag_vs_year_alt.get_yticklabels():
                    tl.set_fontproperties(font)
                    tl.set_fontsize(progress_indicator_fontsize)
                    tl.set_color(progress_tick_color)
                
                for tick in mag_vs_year.yaxis.get_major_ticks():
                    tick.label1.set_alpha(0)
                
                for line in mag_vs_year.xaxis.get_ticklines() + mag_vs_year_alt.yaxis.get_ticklines() + mag_vs_year.yaxis.get_ticklines():
                    line.set_alpha(0)

                #take care of all of the frame line widths
                for spine in mag_vs_year.spines.itervalues():
                    spine.set_lw(progress_frame_width)
                    spine.set_color(progress_frame_color)
                
                mag_vs_year.set_xlabel('magnitude', fontproperties=font, size=progress_indicator_fontsize, color=progress_tick_color)
                mag_vs_year_alt.set_ylabel('year', fontproperties=font, size=progress_indicator_fontsize, color=progress_tick_color)

                #---------------------------------------------------------------
                # Add the progress indicator line
                #---------------------------------------------------------------
                mag_vs_year.axhline(y=progress_indicator_year, lw=progress_indicator_linewidth, c=progress_indicator_line_color)
            
                #---------------------------------------------------------------
                # Add the progress indicator label lines for large events.
                #---------------------------------------------------------------
                label_lines = []
                #if len(event_large_magnitudes) < 10:
                #    label_range = _event_large_magnitude_years
                #else:
                label_range = np.linspace(math.floor(min(_event_years)),math.ceil(max(_event_years)),len(event_large_magnitudes))
                for i, y in enumerate(label_range):
                    m = event_large_magnitudes[i]
                    y1 = _event_large_magnitude_years[i]
                    
                    try:
                        lels = large_event_label_states[event_large_magnitude_evnums[i]]
                    except KeyError:
                        lels = 0

                    the_color = mag_color_map(lels)
                    the_line_color = mag_line_colormap(lels)
                    the_line_width = mag_linewidth + lels * (progress_indicator_linewidth)
                    the_bbox = dict(
                        facecolor=current_mag_face_colormap(lels),
                        linewidth=the_line_width,
                        ec=the_line_color,
                        boxstyle='round4,pad=0.5'
                    )
                    mag_vs_year.text(
                        min(event_magnitudes) - 0.4, y, '{:0.1f}'.format(m),
                        fontproperties=font,
                        fontsize=mag_fontsize,
                        horizontalalignment='right',
                        verticalalignment='center',
                        color=the_color,
                        bbox=the_bbox
                    )
                    
                    label_lines.append(
                        mlines.Line2D(
                            [1.0,-0.01,-0.05,-0.085],
                            [y1/math.ceil(max(_event_years)),
                                y1/math.ceil(max(_event_years)),
                                y/math.ceil(max(_event_years)),
                                y/math.ceil(max(_event_years))
                            ],
                            linewidth=the_line_width,
                            transform=mag_vs_year.transAxes,
                            color=the_line_color,
                            solid_capstyle='round',
                            solid_joinstyle='round'
                        )
                    )
                
                mag_vs_year.lines.extend(label_lines)
                
                #---------------------------------------------------------------
                # Plot the current frequency-magnitude
                #---------------------------------------------------------------
                width_frac = 150.0/pw
                height_frac = 150.0/ph
                left_frac = 100.0/pw
                bottom_frac = 100.0/ph
                
                current_total_events = len(cumulative_magnitudes)
                if current_total_events > 1:
                    cum_freq = {}
                    fm_x = []
                    fm_y = []
                    for num, magnitude in enumerate(sorted(cumulative_magnitudes)):
                        cum_freq['{:0.10f}'.format(magnitude)] = current_total_events - (num + 1)
                    
                    for magnitude in sorted(cum_freq.iterkeys()):
                        fm_x.append(magnitude)
                        fm_y.append(float(cum_freq[magnitude]))
                    mag_vs_freq = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                    mag_vs_freq.semilogy(fm_x, fm_y, color=fm_line_color, linewidth=fm_line_width)
                    mag_vs_freq.set_ylim(bottom=min(fm_y_ticks), top=max(fm_y_ticks))
                    mag_vs_freq.set_xlim(left=round(min(event_magnitudes)), right=round(max(event_magnitudes)))
                    mag_vs_freq.set_yticks(fm_y_ticks)
                    mag_vs_freq.set_xticks(fm_x_ticks)
                    
                    mag_vs_freq_x_ticklabels = ['{:0.1f}'.format(tick) for tick in fm_x_ticks]
                    mag_vs_freq.set_xticklabels(mag_vs_freq_x_ticklabels)
                    
                    for tick in mag_vs_freq.xaxis.get_major_ticks() + mag_vs_freq.yaxis.get_major_ticks():
                        tick.label1.set_fontproperties(font)
                        tick.label1.set_fontsize(fm_fontsize)
                        tick.label1.set_color(fm_label_color)
                        tick.label1.set_alpha(fm_alpha_state)
                    
                    for line in mag_vs_freq.xaxis.get_majorticklines() + mag_vs_freq.yaxis.get_majorticklines() + mag_vs_freq.yaxis.get_minorticklines():
                        line.set_alpha(0)
                    
                    #take care of all of the frame line widths
                    for spine in mag_vs_freq.spines.itervalues():
                        spine.set_lw(fm_frame_width)
                        spine.set_color(fm_frame_color)
                        spine.set_alpha(fm_alpha_state)
                    mag_vs_freq.patch.set_alpha(fm_alpha_state)

                #---------------------------------------------------------------
                # Plot the fault traces.
                #---------------------------------------------------------------
                # Plot the sections
                for sid, sec_trace in fault_traces_latlon.iteritems():
                    trace_Xs, trace_Ys = m4(sec_trace[1], sec_trace[0])
                    
                    try:
                        section_state = section_states[sid]
                    except KeyError:
                        section_state = 0.0

                    m4.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=fault_width + section_state * 2.0, solid_capstyle='round', solid_joinstyle='round')

                #---------------------------------------------------------------
                # Plot the colorbar
                #---------------------------------------------------------------
                left_frac = 70.0/pw
                bottom_frac = (70.0 - cb_height - cb_margin_t)/ph
                width_frac = mw/pw
                height_frac = cb_height/ph
                
                cb_ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
                norm = EFP.norm
                cb = mcolorbar.ColorbarBase(cb_ax, cmap=cmap,
                       norm=norm,
                       orientation='horizontal')
                if field_type == 'displacement':
                    if fringes:
                        cb_title = 'Displacement [m]'
                    else:
                        cb_title = 'Total displacement [m]'

                elif field_type == 'gravity':
                    cb_title = r'Gravity changes [$\mu gal$]'

                cb_ax.set_title(cb_title, fontproperties=font, color=cb_fontcolor, size=cb_fontsize, va='top', ha='left', position=(0,-1.5) )

                for label in cb_ax.xaxis.get_ticklabels():
                    label.set_fontproperties(font)
                    label.set_fontsize(cb_fontsize)
                    label.set_color(cb_fontcolor)
                                      
                for line in cb_ax.xaxis.get_ticklines():
                    line.set_alpha(0)
            
                #---------------------------------------------------------------
                # If the field is a displacement field, draw the look arrows.
                #---------------------------------------------------------------
                if field_type == 'displacement':
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

                    az_arrow_start_x    = 0.5 - (0.8/2.0)*math.sin(EFP.look_azimuth)
                    az_arrow_start_y    = 0.5 - (0.8/2.0)*math.cos(EFP.look_azimuth)
                    az_arrow_dx      = 0.8*math.sin(EFP.look_azimuth)
                    az_arrow_dy      = 0.8*math.cos(EFP.look_azimuth)

                    az_ax.arrow( az_arrow_start_x , az_arrow_start_y, az_arrow_dx, az_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='right', length_includes_head=True, lw=1.0, fc='k' )
                    az_ax.add_line(mlines.Line2D((0.5,0.5), (0.5,0.8), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
                    az_ax.add_patch(mpatches.Arc((0.5,0.5), 0.3, 0.3, theta1=90.0 - EF.convert.rad2deg(EFP.look_azimuth), theta2=90.0, fc='none', lw=1.0, ls='dotted', ec='k'))
                    az_ax.text(1.0, 1.0, 'az = {:0.1f}{}'.format(EF.convert.rad2deg(EFP.look_azimuth),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')

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

                    al_arrow_start_x    = 0.1 + 0.8*math.cos(EFP.look_elevation)
                    al_arrow_start_y    = 0.1 + 0.8*math.sin(EFP.look_elevation)
                    al_arrow_dx      = -0.8*math.cos(EFP.look_elevation)
                    al_arrow_dy      = -0.8*math.sin(EFP.look_elevation)

                    al_ax.arrow( al_arrow_start_x , al_arrow_start_y, al_arrow_dx, al_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='left', length_includes_head=True, lw=1.0, fc='k' )
                    al_ax.add_line(mlines.Line2D((0.1,0.9), (0.1,0.1), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
                    al_ax.add_patch(mpatches.Arc((0.1,0.1), 0.5, 0.5, theta1=0.0, theta2=EF.convert.rad2deg(EFP.look_elevation), fc='none', lw=1.0, ls='dotted', ec='k'))
                    al_ax.text(1.0, 1.0, 'al = {:0.1f}{}'.format(EF.convert.rad2deg(EFP.look_elevation),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')
                #---------------------------------------------------------------
                # If the field is a gravity field change outermost tick labels
                # on colorbar.
                #---------------------------------------------------------------
                else:
                    # Want to change outermost tick labels on colorbar
                    #   from 'VALUE','-VALUE' to '>VALUE' and '<-VALUE'    
                    cb_tick_labs = [item.get_text() for item in cb_ax.get_xticklabels()]
                    cb_tick_labs[0] = '<'+cb_tick_labs[0]
                    cb_tick_labs[-1] = '>'+cb_tick_labs[-1]
                    cb_ax.set_xticklabels(cb_tick_labs)
                
                #---------------------------------------------------------------
                # Save the figure and clear out all matplotlib figures.
                #---------------------------------------------------------------
                fig4.savefig('{}{}.png'.format(frame_images_directory, the_frame), format='png', dpi=plot_resolution)
                
                fig4.clf()
                mplt.close('all')
                gc.collect()

                sys.stdout.write('\033[2K')
                sys.stdout.flush()
            
            #-------------------------------------------------------------------
            # State variables that need to be maintained post-plotting.
            #-------------------------------------------------------------------
            
            # Big magnitude event label fade state
            for k, lels in large_event_label_states.iteritems():
                if lels >= 0.04:
                    large_event_label_states[k] -= 0.04
                else:
                    large_event_label_states[k] = 0.0
            
            # Active section line width
            for sid, section_state in section_states.iteritems():
                if section_state >= 0.04:
                    section_states[sid] -= 0.04
                else:
                    section_states[sid] = 0.0

            # Frequency magnitude plot alpha
            fm_alpha_state += 0.04
            if fm_alpha_state > 1.0:
                fm_alpha_state = 1.0

            sys.stdout.write('\033[1A')
            sys.stdout.write('\033[2K')
            sys.stdout.write('\033[1A\r')
            sys.stdout.write('\033[2K')
        
        #-----------------------------------------------------------------------
        # Create the movie using ffmpeg.
        #-----------------------------------------------------------------------
        
        #ffmpeg -y -r 15 -i f_%d.png -f mp4 -vcodec h264 -pix_fmt yuv420p animation.mp4
        
        proc_args = 'ffmpeg -y -r {fps} -start_number 0 -i {dir}{inc}.png -f mp4 -vcodec h264 -pix_fmt yuv420p {out}animation.mp4'.format(
            fps=int(animation_fps),
            dir=frame_images_directory,
            inc='%d',
            out=output_directory
        )
        proc = subprocess.Popen(proc_args, shell=True)
        proc.wait()




    
#----------------------------------------------------------------------------    
# This method calculates and saves the individual snap shots of the gravity
# fields during the simulation at a specified time before or after a specified
# earthquake. For instance, to compute the characteristic gravity change pattern
# 5 years prior to earthquakes on the southern San Andreas Fault, one needs to
# set buffer=5.0, pre=True, section_filter=(southern SAF sections), and specify
# the event IDs that will be used in this averaging. Then you just pass the 
# plotting objects to the generate_map method to safe the plot.
#---------------------------------------------------------------------------- 

def compute_composite_fields(sim_file, output_directory, event_ids,
    field_type='gravity', fringes=True, padding=0.08, cutoff=None,
    pre=True,buffer=5.0,section_filter=None,
    backslip_only=False,eq_slip_only=False):
    
    # ----------------------------- Initializing --------------------------
    sys.stdout.write('Initializing :: ')
    sys.stdout.flush()
    
    # create the animation dir if needed
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    if not output_directory.endswith('/'):
        output_directory += '/'

    
    # create the animation subdirs if needed
    field_values_directory = '{}field_values/'.format(output_directory)
    frame_images_directory = '{}field_images/'.format(output_directory)
    if not os.path.exists(field_values_directory):
        os.makedirs(field_values_directory)
    if not os.path.exists(frame_images_directory):
        os.makedirs(frame_images_directory)

    #---------------------------------------------------------------------------
    # Open the data file. It needs to stay open while we do the animation.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        events      = VCEvents(sim_data)
        geometry    = VCGeometry(sim_data)
        
        # Get global information about the simulations geometry
        min_lat     = geometry.min_lat
        max_lat     = geometry.max_lat
        min_lon     = geometry.min_lon
        max_lon     = geometry.max_lon
        base_lat    = geometry.base_lat
        base_lon    = geometry.base_lon
        fault_traces= geometry.get_fault_traces()
    
        # ------------------------------------------
        # The blocks in slip_rates define the elements that are used for the plotting
        slip_rates  = geometry.get_slip_rates(section_filter=section_filter,per_year=True)

        # ------------------------------------------
        # Grab rows from event table for each specified event
        avg_event_data  = events.get_event_data_from_evnums(event_ids,['event_magnitude', 'event_year', 'event_number'],section_filter=section_filter)
        
        # Extract relevant data from events inthe time range of the animation
        avg_event_magnitudes    = avg_event_data['event_magnitude']
        avg_event_years         = avg_event_data['event_year']
        avg_event_numbers       = avg_event_data['event_number']
        
        # These are the event years shifted by the buffer length (e.g. 5 years)
        if pre:
            _event_years = [y - buffer for y in avg_event_years]
        else:
            _event_years = [y + buffer for y in avg_event_years]
            
        
        # -------------------------------------------------------------
        # Instantiate the field and the plotter
        # -------------------------------------------------------------
        if field_type == 'displacement':
            EF = vcutils.VCDisplacementField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
            EFP = vcplotutils.VCDisplacementFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)
            EFP.calculate_look_angles(geometry[:])
        elif field_type == 'gravity':
            EF = vcutils.VCGravityField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
            EFP = vcplotutils.VCGravityFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)
        
        
        # -------------------------------------------------------------
        # Loop over each event, and calculate or load the field
        # -------------------------------------------------------------
        for k in range(len(event_ids)):
            evnum               = event_ids[k]
            ev_year            = _event_years[k]
            #event_sections      = geometry.sections_with_elements(event_element_slips[evnum].keys())
            
            
            #-------------------------------------------------------------------    
            # Apply back slip to all elements, up to time of this event
            if not eq_slip_only:
                frame_slips  = {bid:-1.0*ev_year*float(slip_rates[bid]) for bid in slip_rates.keys()}
            else:
                frame_slips  = {bid:0.0 for bid in slip_rates.keys()}

            #--------------------------------------------------------------
            # Set up the elements to evaluate Green's functions
            ele_getter         = itemgetter(*frame_slips.keys())
            frame_element_data = ele_getter(geometry)
            
            #-------------------------------------------------------
            # Grab all events before and including the current event at _ev_year
            if not backslip_only:
                prev_filter         = {'type':'year', 'filter':(0.0,ev_year)}
                all_events          = events.get_event_data(['event_magnitude', 'event_year', 'event_number'], event_range=prev_filter,section_filter=section_filter)

                #-------------------------------------------------------
                # Grab the event element slips for all previous events and add to the back slip            
                for evid in all_events['event_number']:
                    # Grab all event slips up to this event year on specified sections 
                    event_element_slips = events.get_event_element_slips(evid)
                

                    #  Tell John about this step, verify it's ok to leave out these elements
                    for bid in event_element_slips.keys():
                        try:
                            frame_slips[bid] += float(event_element_slips[bid])
                        except KeyError:
                            pass        



            #sys.stdout.write('event {}, year {}\n'.format(evnum, ev_year))

        
            #-------------------------------------------------------------------
            # Calculate all of the data for the current field. 
            #-------------------------------------------------------------------
            if not eq_slip_only and not backslip_only:
                PRE = '{}{}_'.format(field_values_directory, evnum)
            elif eq_slip_only:
                PRE = '{}{}_eq_'.format(field_values_directory, evnum)
            elif backslip_only:
                PRE = '{}{}_back_'.format(field_values_directory, evnum)
            
            sys.stdout.write('\nprocessing event {}'.format(evnum))
            sys.stdout.flush()
            
            sys.stdout.write(', {} elements :: '.format(len(frame_slips.keys())))
            sys.stdout.flush()
                                   

            EF.calculate_field_values(
                    frame_element_data,
                    frame_slips,
                    cutoff=cutoff,
                    save_file_prefix=PRE)
#-------------------------------------------------------------------------------    




#-------------------------------------------------------------------------------
# Plot backslip
#-------------------------------------------------------------------------------
def plot_backslip(sim_file, duration, section_filter=None, field_type='gravity',cutoff=None,padding=0.08,tag=None,fringes=False):
    
    output_directory       = 'backslip_only/'
    field_values_directory = '{}field_values/'.format(output_directory)
    
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        
    if not os.path.exists(field_values_directory):
        os.makedirs(field_values_directory)
        
    
    
    # ----------------------------- Initializing --------------------------
    sys.stdout.write('Initializing plot :: ')
    sys.stdout.flush()
        
    #---------------------------------------------------------------------------
    # Open the data file.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        geometry    = VCGeometry(sim_data)
        events      = VCEvents(sim_data)
        
        # Get global information about the simulations geometry
        min_lat     = geometry.min_lat
        max_lat     = geometry.max_lat
        min_lon     = geometry.min_lon
        max_lon     = geometry.max_lon
        base_lat    = geometry.base_lat
        base_lon    = geometry.base_lon
        fault_traces= geometry.get_fault_traces()
    
        # ------------------------------------------
        # The blocks in slip_rates define the elements that are used for the plotting
        slip_rates         = geometry.get_slip_rates(section_filter=section_filter,per_year=True)
        # Set up the elements to evaluate Green's functions
        ele_getter         = itemgetter(*slip_rates.keys())
        element_data       = ele_getter(geometry)
        
        # Get event information, filter by section if specified
        event_data = events.get_event_data(['event_magnitude', 'event_year', 'event_number'], section_filter=section_filter)

        
    
    # -------------------------------------------------------------
    # Instantiate the field and the plotter
    # -------------------------------------------------------------
    if field_type == 'displacement':
        EF = vcutils.VCDisplacementField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
        EFP = vcplotutils.VCDisplacementFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)
        EFP.calculate_look_angles(geometry[:])
    elif field_type == 'gravity':
        EF = vcutils.VCGravityField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
        EFP = vcplotutils.VCGravityFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)
    
    
    # Apply backslip
    element_slips  = {bid:-1.0*duration*float(slip_rates[bid]) for bid in slip_rates.keys()}

    PRE = '{}{}_'.format(field_values_directory, int(duration)) 

    # Try and load the fields
    field_values_loaded = EF.load_field_values(PRE)
    if field_values_loaded:
        sys.stdout.write('\nloaded {} years of backslip'.format(int(duration)))
    else:
        # If they havent been saved then we need to calculate them
        sys.stdout.write('\nprocessing {} elements :: '.format(len(element_slips.keys())))
        sys.stdout.flush()
            
        EF.calculate_field_values(
                element_data,
                element_slips,
                cutoff=cutoff,
                save_file_prefix=PRE)

    # Make the plot and save it
    generate_map(EF,EFP,fault_traces,fringes,event_data,output_file,field_type='gravity')


#--------------------------------------------------------------------------
def make_composite_field(sim_file, event_ids,field_dir,
    field_type='gravity', fringes=True, padding=0.08, cutoff=None,
    pre=True,buffer=5.0,section_filter=None,
    backslip_only=False,eq_slip_only=False,tag=None):

    out_dir = 'composite_fields/'
    FACTOR  = float(1.0/len(event_ids))


    if not field_dir.endswith('/'):
        field_dir += '/'    
              

    prefix0  = out_dir+'avg_'
    if pre:
        prefix1 = 'pre_'
    else:
        prefix1 = 'post_'
    prefix1 += str(int(buffer[0]))+'yr_'
    
    if backslip_only:
        prefix1 += 'back_'
    if eq_slip_only:
        prefix1 += 'eq_'

       
    if tag is not None:
        prefix1 += tag+'_'   
       
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    #---------------------------------------------------------------------------
    # Open the data file. It needs to stay open while we do the animation.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        geometry    = VCGeometry(sim_data)
        events      = VCEvents(sim_data)
        
        # Get event information, filter by section if specified
        event_data = events.get_event_data(['event_magnitude', 'event_year', 'event_number'], section_filter=section_filter)
        
        # Get global information about the simulations geometry
        min_lat     = geometry.min_lat
        max_lat     = geometry.max_lat
        min_lon     = geometry.min_lon
        max_lon     = geometry.max_lon
        base_lat    = geometry.base_lat
        base_lon    = geometry.base_lon
        fault_traces= geometry.get_fault_traces()
   

    # -------------------------------------------------------------
    # Instantiate the field and the plotter
    # -------------------------------------------------------------
    if field_type == 'displacement':
        EF = vcutils.VCDisplacementField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
        EFP = vcplotutils.VCDisplacementFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)
        EFP.calculate_look_angles(geometry[:])
    elif field_type == 'gravity':
        EF = vcutils.VCGravityField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
        EFP = vcplotutils.VCGravityFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)


    output_file = prefix0+prefix1+'cbar{}.png'.format(EFP.dmc['cbar_max'])

    # -------------------------------------------------------------
    # Add up all the field values that comprise field #1
    # -------------------------------------------------------------
    for evnum in event_ids:
        PRE = field_dir+'field_values/{}_'.format(evnum)
        field_values_loaded = EF.load_field_values(PRE,factor=FACTOR)
        if not field_values_loaded:
            print "Cannot load field values!"
        else:
            sys.stdout.write('\nloaded event {}'.format(evnum))
    sys.stdout.write('\nComposite fields loaded\n')


  
    #------------------------------------
    sys.stdout.write('\nComposite fields complete')
    # Make the plot and save it
    generate_map(EF,EFP,fault_traces,fringes,event_data,output_file,field_type='gravity')



#--------------------------------------------------------------------------
def diff_composite_fields(sim_file, event_ids,field1dir,field2dir,
    field_type='gravity', fringes=True, padding=0.08, cutoff=None,
    pre=(False,True),buffer=(5.0,5.0),section_filter=None,
    backslip_only=(False,False),eq_slip_only=(False,False),tag=None):

    out_dir = 'field_diff/'
    factor1 = float(1.0/len(event_ids[0]))
    factor2 = float(1.0/len(event_ids[1]))

    if not field1dir.endswith('/'):
        output_directory += '/'    
        
    if not field2dir.endswith('/'):
        output_directory += '/'        

    prefix0  = out_dir+'diff_'
    if pre[0]:
        prefix1 = 'pre_'
    else:
        prefix1 = 'post_'
    prefix1 += str(int(buffer[0]))+'yr_'
    
    if pre[1]:
        prefix2 = 'pre_'
    else:
        prefix2 = 'post_'
    prefix2 += str(int(buffer[1]))+'yr_'
    
    if backslip_only[0]:
        prefix1 += 'back_'
    if eq_slip_only[0]:
        prefix1 += 'eq_'
    if backslip_only[1]:
        prefix2 += 'back_'  
    if eq_slip_only[1]:
        prefix2 += 'eq_' 
       
    if tag is not None:
        prefix2 += tag+'_'   
       
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    #---------------------------------------------------------------------------
    # Open the data file. It needs to stay open while we do the animation.
    #---------------------------------------------------------------------------
    with VCSimData() as sim_data:
        # open the simulation data file
        sim_data.open_file(sim_file)
        
        # instantiate the vc classes passing in an instance of the VCSimData
        # class
        geometry    = VCGeometry(sim_data)
        events      = VCEvents(sim_data)
        
        # Get event information, filter by section if specified
        event_data = events.get_event_data(['event_magnitude', 'event_year', 'event_number'], section_filter=section_filter)
        
        # Get global information about the simulations geometry
        min_lat     = geometry.min_lat
        max_lat     = geometry.max_lat
        min_lon     = geometry.min_lon
        max_lon     = geometry.max_lon
        base_lat    = geometry.base_lat
        base_lon    = geometry.base_lon
        fault_traces= geometry.get_fault_traces()
   

    # -------------------------------------------------------------
    # Instantiate the field and the plotter
    # -------------------------------------------------------------
    if field_type == 'displacement':
        EF = vcutils.VCDisplacementField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
        EFP = vcplotutils.VCDisplacementFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)
        EFP.calculate_look_angles(geometry[:])
    elif field_type == 'gravity':
        EF = vcutils.VCGravityField(min_lat, max_lat, min_lon, max_lon, base_lat, base_lon, padding=padding)
        EFP = vcplotutils.VCGravityFieldPlotter(EF.min_lat, EF.max_lat, EF.min_lon, EF.max_lon)


    output_file = prefix0+prefix1+prefix2+'cbar{}.png'.format(EFP.dmc['cbar_max'])

    # -------------------------------------------------------------
    # Add up all the field values that comprise field #1
    # -------------------------------------------------------------
    for evnum in event_ids[0]:
        PRE = field1dir+'field_values/'+'{}_'.format(evnum)
        field_values_loaded = EF.load_field_values(PRE,factor=factor1)
        if not field_values_loaded:
            print "Cannot load field values!"
        else:
            sys.stdout.write('\nloaded event {}'.format(evnum))
    sys.stdout.write('\nComposite 1 loaded\n')


    # -------------------------------------------------------------
    # Add up all the field values that comprise field #2, subtract
    # -------------------------------------------------------------
    for evnum in event_ids[1]:
        PRE = field2dir+'field_values/'+'{}_'.format(evnum)
        field_values_loaded = EF.load_field_values(PRE,factor=factor2,subtract=True)
        if not field_values_loaded:
            print "Cannot load field values!"
        else:
            sys.stdout.write('\nloaded event {}'.format(evnum))
    sys.stdout.write('\nComposite 2 loaded\n')
  
    #------------------------------------
    sys.stdout.write('\nComposite fields complete')
    # Make the plot and save it
    generate_map(EF,EFP,fault_traces,fringes,event_data,output_file,field_type='gravity')
    
    
    
    
#---------------------------------------------------------------------------
def generate_map(EF,EFP,fault_traces,fringes,event_data,output_file,event_sections,field_type='gravity',hi_res=False):

    # Send field values to be plotted
    EFP.set_field(EF) 

  
    sys.stdout.write('\nmap image : ')
    sys.stdout.flush()

    map_image = EFP.create_field_image(fringes=fringes)

    sys.stdout.write('map overlay : ')
    sys.stdout.flush()
    # Convert the fault traces to lat-lon
    fault_traces_latlon = {}
    for secid in fault_traces.iterkeys():
         fault_traces_latlon[secid] = zip(*[(lambda y: (y.lat(),y.lon()))(EF.convert.convert2LatLon(quakelib.Vec3(x[0], x[1], x[2]))) for x in fault_traces[secid]])

    #---------------------------------------------------------------------------
    # Plot all of the geographic info on top of the displacement map image.
    #---------------------------------------------------------------------------
    
    # Grab all of the plot properties that we will need.
    # properties that are fringes dependent
    if fringes and field_type == 'displacement':
        cmap            = EFP.dmc['cmap_f']
        coastline_color = EFP.dmc['coastline_color_f']
        country_color   = EFP.dmc['country_color_f']
        state_color     = EFP.dmc['state_color_f']
        fault_color     = EFP.dmc['fault_color_f']
        map_tick_color  = EFP.dmc['map_tick_color_f']
        map_frame_color = EFP.dmc['map_frame_color_f']
        grid_color      = EFP.dmc['grid_color_f']
        cb_fontcolor    = EFP.dmc['cb_fontcolor_f']
    else:
        cmap            = EFP.dmc['cmap']
        coastline_color = EFP.dmc['coastline_color']
        country_color   = EFP.dmc['country_color']
        state_color     = EFP.dmc['state_color']
        fault_color     = EFP.dmc['fault_color']
        map_tick_color  = EFP.dmc['map_tick_color']
        map_frame_color = EFP.dmc['map_frame_color']
        grid_color      = EFP.dmc['grid_color']
        cb_fontcolor    = EFP.dmc['cb_fontcolor']
    
    # properties that are not fringes dependent
    boundary_width  = EFP.dmc['boundary_width']
    coastline_width = EFP.dmc['coastline_width']
    country_width   = EFP.dmc['country_width']
    state_width     = EFP.dmc['state_width']
    river_width     = EFP.dmc['river_width']
    fault_width     = EFP.dmc['fault_width']
    map_frame_width = EFP.dmc['map_frame_width']
    map_fontsize    = EFP.dmc['map_fontsize']
    arrow_inset     = EFP.dmc['arrow_inset']
    arrow_fontsize  = EFP.dmc['arrow_fontsize']
    cb_fontsize     = EFP.dmc['cb_fontsize']
    cb_height       = EFP.dmc['cb_height']
    cb_margin_t     = EFP.dmc['cb_margin_t']
    grid_width      = EFP.dmc['grid_width']
    num_grid_lines  = EFP.dmc['num_grid_lines']
    font            = EFP.dmc['font']
    font_bold       = EFP.dmc['font_bold']

    map_resolution  = EFP.dmc['map_resolution']
    map_projection  = EFP.dmc['map_projection']
    plot_resolution = EFP.dmc['plot_resolution']

    # The sizing for the image is tricky. The aspect ratio of the plot is fixed,
    # so we cant set all of margins to whatever we want. We will set the anchor
    # to the top, left margin position. Then scale the image based on the
    # bottom/right margin, whichever is bigger.
    
    mw = EF.lons_1d.size
    mh = EF.lats_1d.size

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

    if hi_res:
        fig_res = plot_resolution*4.0
    else:
        fig_res = plot_resolution

    fig4 = mplt.figure(figsize=(pwi, phi), dpi=fig_res)

    #---------------------------------------------------------------------------
    # m4, fig4 is all of the boundary data.
    #---------------------------------------------------------------------------
    m4 = Basemap(
        llcrnrlon=EF.min_lon,
        llcrnrlat=EF.min_lat,
        urcrnrlon=EF.max_lon,
        urcrnrlat=EF.max_lat,
        lat_0=(EF.max_lat+EF.min_lat)/2.0,
        lon_0=(EF.max_lon+EF.min_lon)/2.0,
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
    parallels = np.linspace(EFP.lats_1d.min(), EFP.lats_1d.max(), num_grid_lines+1)
    m4_parallels = m4.drawparallels(parallels, labels=[1,0,0,0], color=grid_color, fontproperties=font, fmt='%.2f', linewidth=grid_width, dashes=[1, 10], fontsize=map_fontsize)
    
    # draw meridians
    meridians = np.linspace(EFP.lons_1d.min(), EFP.lons_1d.max(), num_grid_lines+1)
    m4_meridians = m4.drawmeridians(meridians, labels=[0,0,1,0], color=grid_color, fontproperties=font, fmt='%.2f', linewidth=grid_width, dashes=[1, 10], fontsize=map_fontsize)

    if field_type == 'displacement':
        box_size = 70.0
        # draw the azimuth look arrow
        az_width_frac    = box_size/pw
        az_height_frac   = box_size/ph
        az_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
        az_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac)/ph
        az_ax = fig4.add_axes((az_left_frac,az_bottom_frac,az_width_frac,az_height_frac))

        az_ax.set_xlim((0,1.0))
        az_ax.set_ylim((0,1.0))
        for item in az_ax.yaxis.get_ticklabels() + az_ax.xaxis.get_ticklabels() + az_ax.yaxis.get_ticklines() + az_ax.xaxis.get_ticklines():
            item.set_alpha(0)

        az_arrow_start_x    = 0.5 - (0.8/2.0)*math.sin(EFP.look_azimuth)
        az_arrow_start_y    = 0.5 - (0.8/2.0)*math.cos(EFP.look_azimuth)
        az_arrow_dx      = 0.8*math.sin(EFP.look_azimuth)
        az_arrow_dy      = 0.8*math.cos(EFP.look_azimuth)

        az_ax.arrow( az_arrow_start_x , az_arrow_start_y, az_arrow_dx, az_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='right', length_includes_head=True, lw=1.0, fc='k' )
        az_ax.add_line(mlines.Line2D((0.5,0.5), (0.5,0.8), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
        az_ax.add_patch(mpatches.Arc((0.5,0.5), 0.3, 0.3, theta1=90.0 - EF.convert.rad2deg(EFP.look_azimuth), theta2=90.0, fc='none', lw=1.0, ls='dotted', ec='k'))
        az_ax.text(1.0, 1.0, 'az = {:0.1f}{}'.format(EF.convert.rad2deg(EFP.look_azimuth),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')

        # draw the altitude look arrow
        al_width_frac    = box_size/pw
        al_height_frac   = box_size/ph
        al_left_frac     = (70.0 + mw - arrow_inset - pw*az_width_frac)/pw
        al_bottom_frac   = (70.0 + mh - arrow_inset - ph*az_height_frac - ph*al_height_frac)/ph
        al_ax = fig4.add_axes((al_left_frac,al_bottom_frac,al_width_frac,al_height_frac))

        al_ax.set_xlim((0,1.0))
        al_ax.set_ylim((0,1.0))
        for item in al_ax.yaxis.get_ticklabels() + al_ax.xaxis.get_ticklabels() + al_ax.yaxis.get_ticklines() + al_ax.xaxis.get_ticklines():
            item.set_alpha(0)

        al_arrow_start_x    = 0.1 + 0.8*math.cos(EFP.look_elevation)
        al_arrow_start_y    = 0.1 + 0.8*math.sin(EFP.look_elevation)
        al_arrow_dx      = -0.8*math.cos(EFP.look_elevation)
        al_arrow_dy      = -0.8*math.sin(EFP.look_elevation)

        al_ax.arrow( al_arrow_start_x , al_arrow_start_y, al_arrow_dx, al_arrow_dy, head_width=0.1, head_length= 0.1, overhang=0.1, shape='left', length_includes_head=True, lw=1.0, fc='k' )
        al_ax.add_line(mlines.Line2D((0.1,0.9), (0.1,0.1), lw=1.0, ls=':', c='k', dashes=(2.0,1.0)))
        al_ax.add_patch(mpatches.Arc((0.1,0.1), 0.5, 0.5, theta1=0.0, theta2=EF.convert.rad2deg(EFP.look_elevation), fc='none', lw=1.0, ls='dotted', ec='k'))
        al_ax.text(1.0, 1.0, 'al = {:0.1f}{}'.format(EF.convert.rad2deg(EFP.look_elevation),r'$^{\circ}$'), fontproperties=font_bold, size=arrow_fontsize, ha='right', va='top')
        
        # draw the box with the magnitude
        mag_width_frac    = box_size/pw
        mag_height_frac   = 15.0/ph    # originally 10.0/ph
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
    for sid, sec_trace in fault_traces_latlon.iteritems():
        trace_Xs, trace_Ys = m4(sec_trace[1], sec_trace[0])
        
        if sid in event_sections:
            linewidth = fault_width + 2.5
        else:
            linewidth = fault_width

        m4.plot(trace_Xs, trace_Ys, color=fault_color, linewidth=linewidth, solid_capstyle='round', solid_joinstyle='round')

    #plot the cb
    left_frac = 70.0/pw
    bottom_frac = (70.0 - cb_height - cb_margin_t)/ph
    width_frac = mw/pw
    height_frac = cb_height/ph
    
    cb_ax = fig4.add_axes((left_frac,bottom_frac,width_frac,height_frac))
    norm = EFP.norm
    cb = mcolorbar.ColorbarBase(cb_ax, cmap=cmap,
           norm=norm,
           orientation='horizontal')
    if field_type == 'displacement':
        if fringes:
            cb_title = 'Displacement [m]'
        else:
            cb_title = 'Total displacement [m]'

    elif field_type == 'gravity' or field_type == 'dilat_gravity':
        cb_title        = r'Gravity changes [$\mu gal$]'

        # If contour plot, match colorbar ticks to contour values
        if EFP.contours is not None:
            cb.set_ticks(EFP.contours)

        # Make first and last ticks on colorbar be <MIN and >MAX
        cb_tick_labs    = [item.get_text() for item in cb_ax.get_xticklabels()]
        cb_tick_labs[0] = '<'+cb_tick_labs[0]
        cb_tick_labs[-1]= '>'+cb_tick_labs[-1]
        cb_ax.set_xticklabels(cb_tick_labs)

    cb_ax.set_title(cb_title, fontproperties=font, color=cb_fontcolor, size=cb_fontsize, va='top', ha='left', position=(0,-1.5) )

    for label in cb_ax.xaxis.get_ticklabels():
        label.set_fontproperties(font)
        label.set_fontsize(cb_fontsize)
        label.set_color(cb_fontcolor)
    for line in cb_ax.xaxis.get_ticklines():
        line.set_alpha(0)


    fig4.savefig(output_file, format='png', dpi=fig_res)

    sys.stdout.write('\nPlot saved: {}'.format(output_file))
    sys.stdout.write('\ndone\n')
    sys.stdout.flush()
#---------------------------------------------------------------------------
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




