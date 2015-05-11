#!/usr/bin/env python
import matplotlib.pyplot as mplt
import matplotlib.font_manager as mfont
import matplotlib.colors as mcolor
import matplotlib.colorbar as mcolorbar
import matplotlib.lines as mlines
from mpl_toolkits.basemap import Basemap
from PIL import Image

import numpy as np

import math
import gc
import multiprocessing

import quakelib

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
# standard plotting routine for scaling plots
#-------------------------------------------------------------------------------
def standard_plot(output_file, x, y, legend_loc='best', scale=False, **kwargs):
    add_lines = kwargs.get('add_lines')
    axis_labels = kwargs.get('axis_labels')
    plot_label = kwargs.get('plot_label')
    connect_points = kwargs.get('connect_points')
    axis_format = kwargs.get('axis_format')
    
    freq_mag = False
    
    if output_file is not None:
        plot_format = output_file.split('.')[-1]
    else:
        plot_format = 'png'

    if plot_format != 'png' and plot_format != 'pdf' and plot_format != 'dat' and plot_format != 'tsv':
            raise vcexceptions.PlotFormatNotSupported(plot_format)
    elif plot_format == 'png' or plot_format == 'pdf':
    #---------------------------------------------------------------------------
    # plot the data using matplotlib
    #---------------------------------------------------------------------------
        
        #-----------------------------------------------------------------------
        # set up plot dimensions (all values in pixels)
        #-----------------------------------------------------------------------
        # the full image width
        imw = 501.0+15
        # the full image height
        imh = 501.0 
        # the left margin and bottom margin
        if axis_labels is not None:
            lm = 45.0+15
            bm = 45.0
        else:
            lm = 10.0
            bm = 10.0
        # the right margin
        rm = 10.0 
        # the top margin
        if plot_label is not None:
            tm = 30.0
        else:
            tm = 10.0
        # the plot resolution
        res = 72.0

        #-----------------------------------------------------------------------
        # set up the plot styles
        #-----------------------------------------------------------------------
        # the marker face color for the main plot
        mfc_main = (0,0,0,1)
        # the marker size for the main plot
        ms_main = 3.0
        # the marker type for the main plot
        mt_main = 'o'
        # the line style for the main plot
        if connect_points is not None:
            ls_main = '--'
        else:
            ls_main = 'None'
        # the line color for the main plot
        c_main = (0,0,0,1)
        # styles for additional sub plots
        if add_lines is not None:
            # the line style for the extra lines
            ls_extra = '-'
            # the line weight for the extra lines
            lw_extra = 4.0
            # the line color for the extra lines
            c_extra = (0.5,0.5,0.5,1)
    
        #-----------------------------------------------------------------------
        # set up the plot fonts
        #-----------------------------------------------------------------------
        ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=15)
        framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=16)
        legendfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=16)
        if len(plot_label) > 100:
            title_size = 9
        elif len(plot_label) <= 100 and len(plot_label) > 89:
            title_size = 10
        elif len(plot_label) <= 89 and len(plot_label) > 50:
            title_size = 11
        else:
            title_size = 12
        titlefont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=title_size)
        
        #-----------------------------------------------------------------------
        # do the plot
        #-----------------------------------------------------------------------
        # calculate the final dimensions and create the figure and axis
        imwi = imw/res
        imhi = imh/res
        pw = imw - lm - rm
        ph = imh - tm - bm
        fig = mplt.figure(figsize=(imwi, imhi), dpi=res)
        the_ax = fig.add_axes((lm/imw, bm/imh, pw/imw, ph/imh))
        
        # the main plotting routines
        eval('the_ax.{}(x, y, ls=ls_main, mfc=mfc_main, ms=ms_main, marker=mt_main, c=c_main)'.format(axis_format), locals())
        # plot any additional lines
        if add_lines is not None:
            for line in add_lines:
                # Catch if it's a freq-mag plot, need to adjust plotting bounds
                if line['label'] == 'UCERF2':
                    freq_mag = True
                    
                try:
                    ls_extra = line['ls']
                except KeyError:
                    pass
                try:
                    lw_extra = line['lw']
                except KeyError:
                    pass
                try:
                    c_extra = line['c']
                except KeyError:
                    pass
                try:
                    yerr=line['y_error']
                    add_axis_format = 'errorbar'
                    y_error_key = 'yerr=line[\'y_error\'],'
                except KeyError:
                    y_error_key = ''
                    add_axis_format = axis_format
                
                eval('the_ax.{0}(line[\'x\'], line[\'y\'], ls=ls_extra, {1} lw=lw_extra, c=c_extra, label=line[\'label\'])'.format(add_axis_format, y_error_key), locals())
        
        # set the fonts for the tick labels
        for label in the_ax.xaxis.get_ticklabels()+the_ax.yaxis.get_ticklabels():
            label.set_fontproperties(ticklabelfont)
        
        # label the axes
        if axis_labels is not None:
            the_ax.set_ylabel(axis_labels['y'], fontproperties=framelabelfont)
            the_ax.set_xlabel(axis_labels['x'], fontproperties=framelabelfont)

        # label the plot
        if plot_label is not None:
            the_ax.set_title(plot_label, fontproperties=titlefont, x=0, y=1, ha='left')

        # clean up the final plot
        the_ax.autoscale_view(tight=True)
        
        # create a legend if we have extra lines
        if add_lines is not None:
            the_ax.legend(prop=legendfont, loc=legend_loc)
            
        # Schultz: adding this parameter to zoom into only the x/y data range
        if scale:
            if freq_mag:
                the_ax.set_xlim(min(x)*.99,max(x))
                the_ax.set_ylim(min(y),max(y)*5)
            else:
                the_ax.set_xlim(min(x),max(x))
                the_ax.set_ylim(min(y),max(y))
            
        
        if output_file is not None:
            # save the plot
            mplt.savefig(output_file, format=plot_format, dpi=res)
    elif plot_format == 'dat':
        f = open(output_file,'w')

        for i, xi in enumerate(x):
            f.write('{mag} {freq}\r'.format(mag=xi, freq=y[i]))

        f.close()
    elif plot_format == 'tsv':
        f = open(output_file,'w')

        for i, xi in enumerate(x):
            f.write('{mag}\t{freq}\r'.format(mag=xi, freq=y[i]))

        f.close()

#-------------------------------------------------------------------------------
# a class to manage the space-time plot
#-------------------------------------------------------------------------------
class VCSpaceTimePlot(object):
    #---------------------------------------------------------------------------
    # y_axis_data_size is the number of years in the plot, or if event_time is
    # true it is the number of events.
    #---------------------------------------------------------------------------
    def __init__(self, output_file, x_axis_data_size, y_axis_data_size, max_depth, min_mag, max_mag, start_year, max_label_len, event_time=False):
        # store the output file
        self.output_file = output_file
        # store the max depth. this is in number of blocks
        self.max_depth = max_depth
        # store the x and y axis data size
        self.x_axis_data_size = x_axis_data_size
        self.y_axis_data_size = y_axis_data_size
        # store the min and max magnitudes
        self.max_mag = max_mag
        self.min_mag = min_mag
        # calculate the slope of the magnitude color ramp
        self.mag_slope = 1.0/(self.max_mag - self.min_mag)
        # store the max label length used to calculate the top margin
        self.max_label_len = max_label_len
        # set the start and end year
        self.start_year = start_year
        self.end_year = self.start_year + self.y_axis_data_size

        #-----------------------------------------------------------------------
        # set up the plot fonts
        #-----------------------------------------------------------------------
        self.ticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
        self.sectionlabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
        self.framelabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=12)
        self.legendtitlefont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=10)
        self.legendticklabelfont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=9)
        self.titlefont = mfont.FontProperties(family='Arial', style='normal', variant='normal', size=14)

        #-----------------------------------------------------------------------
        # set up plot dimensions (all values in pixels)
        #-----------------------------------------------------------------------
        # the event line width and pixels per "y" are based on the y_axis_data_size
        self.elw = 1.0
        self.ppy = 1.0
        if math.ceil(self.y_axis_data_size) < 512 and math.ceil(self.y_axis_data_size) >= 256:
            if event_time:
                self.elw = 2.0
            else:
                self.ppy = 2.0
        elif math.ceil(self.y_axis_data_size) < 256:
            if event_time:
                self.elw = 4.0
            else:
                self.ppy = 4.0
        
        # the event line width and pixels per "x" are based on the x_axis_data_size
        self.ppx = 1.0
        if math.ceil(self.x_axis_data_size) < 640 and math.ceil(self.x_axis_data_size) >= 320:
            self.ppx = 2.0
        elif math.ceil(self.x_axis_data_size) < 320:
            self.ppx = 4.0
        
        print math.ceil(self.x_axis_data_size)
        # the margins
        self.lm = 50.0
        self.bm = 10.0
        self.rm = 10.0
        
        # the colorbar dimensions
        self.cbw = 200.0
        self.cbh = 10.0
        
        # the padding between the section names and the top of the plot frame
        self.snp = 40.0
        
        # the plot resolution
        self.res = 72.0
        
        #-----------------------------------------------------------------------
        # plot dimensions that are calculated based on the settings above
        #-----------------------------------------------------------------------
        # if we are using event time the pixels per "y" is just the event line
        # width
        if event_time:
            self.ppy = self.elw

        # the height of the plot area depends on the pixels per "y" value
        self.ph = self.ppy * math.ceil(self.y_axis_data_size)

        # the width of the plot area depends on the pixels per "x" value
        self.pw = self.ppx * math.ceil(self.x_axis_data_size)

        # the top margin is based on the size of the section labels and their
        # placement
        self.tm = self.snp + self.sectionlabelfont.get_size() * 0.75 * self.max_label_len + self.cbh + 20

        # the total image height
        self.imh = self.ph + self.tm + self.bm
        # the total image width
        self.imw = self.pw + self.lm + self.rm
        
        #-----------------------------------------------------------------------
        # set up the plot styles
        #-----------------------------------------------------------------------
        # the color map for the plot
        self.cmap = mplt.get_cmap('autumn_r',lut=1000)
        # the color of the vertical section demarcation lines
        self.sdlc = '0.5'
        # the alpha of the vertical section demarcation lines
        self.sdla = 0.5
        # the width of the vertical section demarcation lines
        self.sdlw = 0.5
        # the color of the section label lines
        self.sllc = '0.5'
        # the alpha of the section label lines
        self.slla = 0.5
        # the width of the section label lines
        self.sllw = 0.5

        #-----------------------------------------------------------------------
        # start the plot
        #-----------------------------------------------------------------------
        # calculate the final dimensions and create the figure and axis
        self.imwi = self.imw/self.res
        self.imhi = self.imh/self.res
        self.fig = mplt.figure(figsize=(self.imwi, self.imhi), dpi=self.res)
        self.the_ax = self.fig.add_axes((self.lm/self.imw, self.bm/self.imh, self.pw/self.imw, self.ph/self.imh))
    
    #---------------------------------------------------------------------------
    # Adds event to the plot by plotting them directly instead of storing them.
    # This behavior may have to change to storing when running in parallel.
    #---------------------------------------------------------------------------
    def add_event(self, event_number, event_year, event_magnitude, event_line):
        # initilize state variables
        last_value = None
        last_slot = -1
        in_line = False
        # the color of the event line based on the magnitude
        (r,g,b,a) = self.cmap(self.mag_slope * (event_magnitude - self.min_mag))
        # The slots are the indicies of the non-zero line elements from the
        # event line array. These correspond to x data locations where there is
        # a line.
        slots = np.nonzero(event_line)[0]
        # Go through the slots drawing lines where appropriate.
        for k, slot in enumerate(slots):
            value = event_line[slot]
            if (slot - 1) != last_slot and k != 0:
                # we have jumped a gap
                # draw the last line
                self.the_ax.axhline(
                    y=event_year,
                    xmin=float(line_start)/float(self.x_axis_data_size),
                    xmax=float(line_end + 1)/float(self.x_axis_data_size),
                    #alpha=float(last_value)/float(self.max_depth),
                    linewidth=self.elw,
                    color=(r,g,b)
                )
                in_line = False
                # start a new line
                in_line = True
                line_start = slot
                line_end = slot
            else:
                # we are in a contiguous group
                if last_value != event_line[slot]:
                    if in_line:
                        # draw the last line
                        self.the_ax.axhline(
                            y=event_year,
                            xmin=float(line_start)/float(self.x_axis_data_size),
                            xmax=float(line_end + 1)/float(self.x_axis_data_size),
                            #alpha=float(last_value)/float(self.max_depth),
                            linewidth=self.elw,
                            color=(r,g,b)
                        )
                        in_line = False
                    # start line
                    in_line = True
                    line_start = slot
                    line_end = slot
                else:
                    # add to current line
                    line_end = slot
            
            if k == len(slots) - 1:
                # draw the last line
                self.the_ax.axhline(
                    y=event_year,
                    xmin=float(line_start)/float(self.x_axis_data_size),
                    xmax=float(line_end + 1)/float(self.x_axis_data_size),
                    #alpha=float(value)/float(self.max_depth),
                    linewidth=self.elw,
                    color=(r,g,b)
                )
                in_line = False
            last_value = value
            last_slot = slot
    
    #---------------------------------------------------------------------------
    # Add the section labels and the section demarcation lines.
    #---------------------------------------------------------------------------
    def add_section_labels(self, section_offsets, section_info):
        # an array to store all of the label lines
        label_lines = []
        
        # The sorted section ids. Sections are printed in ascending order (by
        # section id) from left to right.
        section_ids = sorted(section_offsets.keys())
        
        # Go through each section plotting the demarcation line, the section
        # label, and the section label line.
        for snum, sid in enumerate(section_ids):
            # Plot the demarcation line (except at the zero offset)
            if section_offsets[sid] != 0:
                self.the_ax.axvline(
                    x=section_offsets[sid],
                    linewidth=self.sdlw,
                    color=self.sdlc,
                    alpha=self.sdla
                )
            
            # Set the text location. This will evenly distribute the text across
            # the length of the x-axis.
            x_text_loc = (float(snum)+0.5)/float(len(section_ids))
            
            # Find the horizontal center of each section region
            if snum == len(section_ids) - 1:
                x_fault_loc = (float(section_offsets[sid] + self.x_axis_data_size)/2.0)/float(self.x_axis_data_size)
            else:
                x_fault_loc = (float(section_offsets[sid] + section_offsets[section_ids[snum + 1]])/2.0)/float(self.x_axis_data_size)
            
            # Plot the section label
            self.the_ax.text(
                x_text_loc, (self.ph + self.snp)/self.ph,
                '{} {}'.format(sid,section_info[sid]['name']),
                horizontalalignment='center',
                verticalalignment='bottom',
                rotation=90,
                fontproperties=self.sectionlabelfont,
                transform=self.the_ax.transAxes
            )

            # Create the section line
            line_xs = [x_text_loc, x_text_loc, x_fault_loc, x_fault_loc]
            line_ys = [
                (self.ph + self.snp - 5)/self.ph,
                (self.ph + self.snp * ( 3.0/4.0 ))/self.ph,
                (self.ph + self.snp * ( 1.0/4.0 ))/self.ph,
                1
            ]
            label_lines.append(mlines.Line2D(
                line_xs, line_ys,transform=self.the_ax.transAxes, solid_capstyle='round',
                    alpha=self.slla,
                    linewidth=self.sllw,
                    color=self.sllc
            ))
        
        # Add the section label lines to the plot.
        self.the_ax.lines.extend(label_lines)
    
    #---------------------------------------------------------------------------
    # Add the title.
    #---------------------------------------------------------------------------
    def add_title(self, title):
        self.the_ax.set_title(
            title,
            position=((10-self.lm)/self.pw,(self.ph + self.tm - 27)/self.ph),
            ha='left',
            fontproperties=self.titlefont,
        )
    
    #---------------------------------------------------------------------------
    # Finish off the plot and save it.
    #---------------------------------------------------------------------------
    def plot(self):
        # Set the axis limits
        self.the_ax.set_xlim((0, self.x_axis_data_size))
        self.the_ax.set_ylim((self.start_year, self.end_year))
        
        # Create the colorbar
        cb_ax    = self.fig.add_axes([(self.lm + self.pw - self.cbw)/self.imw, (self.imh - self.cbh - 10)/self.imh, self.cbw/self.imw, self.cbh/self.imh])
        norm = mcolor.Normalize(vmin=self.min_mag, vmax=self.max_mag)
        cb = mcolorbar.ColorbarBase(cb_ax, cmap=self.cmap,
                   norm=norm,
                   orientation='horizontal')
        ticks = map(float, [self.min_mag,(self.min_mag+ self.max_mag)/2.0, self.max_mag])
        cb.set_ticks([ticks[0], ticks[1], ticks[2]])
        cb.set_ticklabels(['%.2f'%ticks[0],'%.2f'%ticks[1], '%.2f'%ticks[2]])
        
        # Style and cleanup the colorbar ticks
        for label in cb_ax.xaxis.get_ticklabels():
            label.set_fontproperties(self.legendticklabelfont)
        for line in cb_ax.xaxis.get_ticklines():
            line.set_alpha(0)
        
        # Set the colorbar label
        cb.ax.set_xlabel('Magnitude',position=(1,0), ha='right', fontproperties=self.legendtitlefont)
        
        # Style and cleanup the plot ticks
        for label in self.the_ax.yaxis.get_ticklabels():
            label.set_fontproperties(self.ticklabelfont)
        for label in self.the_ax.xaxis.get_ticklabels():
            label.set_alpha(0)
        for line in self.the_ax.xaxis.get_ticklines():
            line.set_alpha(0)
        
        # Set the plot y-axis label
        self.the_ax.set_ylabel('Year',fontproperties=self.framelabelfont)
        
        if self.output_file is not None:
            # Get the plot format and save the file
            plot_format = self.output_file.split('.')[-1]
            if plot_format != 'png' and plot_format != 'pdf':
                raise vcexceptions.PlotFormatNotSupported(plot_format)
            else:
                self.fig.savefig(self.output_file, format=plot_format, dpi=self.res)

#-------------------------------------------------------------------------------
# A class for plotting spacetime events in parallel. This is currently under
# development.
#-------------------------------------------------------------------------------
# TODO: Figure out a way to plot in parallel.
class SpaceTimePlotter(multiprocessing.Process):
    def __init__(self, space_time_plot_params, work_queue, result_queue):
 
        # job management stuff
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        
        #the path to the sim file
        self.stp_params = space_time_plot_params
        
        super(SpaceTimePlotter, self).__init__()
    
    def run(self):
    
        stp = VCSpaceTimePlot(
            self.stp_params['output_file'],
            self.stp_params['x_axis_data_size'],
            self.stp_params['y_axis_data_size'],
            self.stp_params['max_depth'],
            self.stp_params['min_mag'],
            self.stp_params['max_mag'],
            self.stp_params['start_year'],
            self.stp_params['max_label_len'],
            lite_init=True
        )
        
        geometry = self.stp_params['geometry']
        x_data_size = self.stp_params['x_axis_data_size']
        min_depth = self.stp_params['max_depth']
        section_offsets = self.stp_params['section_offsets']
        
        while not self.kill_received:
            # get a task
            try:
                event_data = self.work_queue.get_nowait()
            except Queue.Empty:
                break
        
            # do the processing
            print 'processing events {} - {}'.format(event_data['event_number'][0],event_data['event_number'][-1])
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
            self.result_queue.put(stp)
        
        #self.sim_file.close()

#-------------------------------------------------------------------------------
# A class to handle plotting event displacement fields
#-------------------------------------------------------------------------------
class VCDisplacementFieldPlotter(object):
    def __init__(self, min_lat, max_lat, min_lon, max_lon, map_res='i', map_proj='cyl'):
        self.look_azimuth = None
        self.look_elevation = None
        self.wavelength = 0.03
        
        self.norm = None
        
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
            'map_fontsize':         36.0, #12
            'arrow_inset':          18.0, #10
            'arrow_fontsize':       14.0, # 9
            'cb_fontsize':          20.0, #10
            'cb_fontcolor':         '#000000',
            'cb_fontcolor_f':       '#000000',
            'cb_height':            20.0,
            'cb_margin_t':          4.0
        }
        
        #-----------------------------------------------------------------------
        # m1, fig1 is the oceans and the continents. This will lie behind the
        # masked data image.
        #-----------------------------------------------------------------------
        self.m1 = Basemap(
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
        #-----------------------------------------------------------------------
        # m2, fig2 is the plotted deformation data.
        #-----------------------------------------------------------------------
        self.m2 = Basemap(
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
        #-----------------------------------------------------------------------
        # m3, fig3 is the ocean land mask.
        #-----------------------------------------------------------------------
        self.m3 = Basemap(
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
            
    def set_field(self, field):
        self.lons_1d = field.lons_1d
        self.lats_1d = field.lats_1d
        self.dX = field.dX
        self.dY = field.dY
        self.dZ = field.dZ


    #---------------------------------------------------------------------------
    # Returns a PIL image of the masked displacement map using the current
    # values of the displacements. This map can then be combined into a still
    # or used as part of an animation.
    #---------------------------------------------------------------------------
    def create_field_image(self, fringes=True):
        
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
        
        #print mw, mh
        
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
        
        
        # make sure the values are located at the correct location on the map
        dMags_transformed = self.m2.transform_scalar(dMags, self.lons_1d, self.lats_1d, self.lons_1d.size, self.lats_1d.size)
        
        # prepare the colors for the plot and do the plot
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
            if self.norm is None:
                self.norm = mcolor.Normalize(vmin=0, vmax=self.wavelength)
            im = self.m2.imshow(dMags_colors, interpolation='spline36')
        else:
            dMags_colors = np.empty(dMags_transformed.shape)
            non_zeros = dMags_transformed.nonzero()
            #vmin = np.amin(np.fabs(dMags_transformed[non_zeros]))
            dMags_colors.fill(5e-4)
            dMags_colors[non_zeros] = np.fabs(dMags_transformed[non_zeros])
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
            if self.norm is None:
                self.norm = mcolor.LogNorm(vmin=5e-4, vmax=mod_vmax, clip=True)
            im = self.m2.imshow(dMags_colors, cmap=cmap, norm=self.norm)
        
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
        
        # Clear all three figures
        fig1.clf()
        fig2.clf()
        fig3.clf()
        mplt.close('all')
        gc.collect()
        
        # The final composited image.
        return  Image.composite(im1, im2, mask)

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
        
    # --------------------------------------------------------------------------    
    def set_look_angles(self, azimuth, elevation):
        self.look_azimuth = azimuth
        self.look_elevation = elevation

        
        
        

#-------------------------------------------------------------------------------
# A class to handle plotting event gravity fields
#-------------------------------------------------------------------------------
class VCGravityFieldPlotter(object):
    def __init__(self, min_lat, max_lat, min_lon, max_lon, map_res='i', map_proj='cyl', contours=None):
        
        self.norm = None
        self.contours = contours
        
        #-----------------------------------------------------------------------
        # Gravity map configuration
        #-----------------------------------------------------------------------
        self.dmc = {
            'font':               mfont.FontProperties(family='Arial', style='normal', variant='normal', weight='normal'),
            'font_bold':          mfont.FontProperties(family='Arial', style='normal', variant='normal', weight='bold'),
            'cmap':               mplt.get_cmap('seismic'),
        #water
            'water_color':          '#4eacf4',
        #map boundaries
            'boundary_color':       '#000000',
            'boundary_width':       1.0,
            'coastline_color':      '#000000',
            'coastline_width':      1.0,
            'country_color':        '#000000',
            'country_width':        1.0,
            'state_color':          '#000000',
            'state_width':          1.0,
        #rivers
            'river_width':          0.25,
        #faults
            'fault_color':          '#000000',
            'event_fault_color':    '#ff0000',
            'fault_width':          0.5,
        #lat lon grid
            'grid_color':           '#000000',
            'grid_width':           0.0,
            'num_grid_lines':       5,
        #map props
            'map_resolution':       map_res,
            'map_projection':       map_proj,
            'plot_resolution':      72.0,
            'map_tick_color':       '#000000',
            'map_frame_color':      '#000000',
            'map_frame_width':      1,
        #map_fontsize = 12
            'map_fontsize':         16.0,   # 12   THIS IS BROKEN
            'arrow_inset':          14.0, # 10
            'arrow_fontsize':       12.0, #  9
        #cb_fontsize = 12
            'cb_fontsize':          20.0, # 12
            'cb_fontcolor':         '#000000',
            'cb_height':            20.0,
            'cb_margin_t':          2.0, # 10
         #min/max gravity change labels for colorbar (in microgals)
            'cbar_min':             -20,
            'cbar_max':             20
        }
        
        #-----------------------------------------------------------------------
        # m1, fig1 is the oceans and the continents. This will lie behind the
        # masked data image.
        #-----------------------------------------------------------------------
        self.m1 = Basemap(
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
        #-----------------------------------------------------------------------
        # m2, fig2 is the plotted deformation data.
        #-----------------------------------------------------------------------
        self.m2 = Basemap(
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
        #-----------------------------------------------------------------------
        # m3, fig3 is the ocean land mask.
        #-----------------------------------------------------------------------
        self.m3 = Basemap(
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

    def set_field(self, field):
        self.lons_1d = field.lons_1d
        self.lats_1d = field.lats_1d
        self.dG = field.dG


    #---------------------------------------------------------------------------
    # Returns a PIL image of the masked displacement map using the current
    # values of the displacements. This map can then be combined into a still
    # or used as part of an animation.
    #---------------------------------------------------------------------------
    def create_field_image(self, fringes=True):
        
        #-----------------------------------------------------------------------
        # Set all of the plotting properties
        #-----------------------------------------------------------------------
    
        cmap            = self.dmc['cmap']
        water_color     = self.dmc['water_color']
        boundary_color  = self.dmc['boundary_color']
        land_color      = cmap(0.5)
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
        
        # make sure the values are located at the correct location on the map
        dG_transformed = self.m2.transform_scalar(self.dG, self.lons_1d, self.lats_1d, self.lons_1d.size, self.lats_1d.size)
        
        if self.norm is None:
            #self.norm = mcolor.Normalize(vmin=np.amin(dG_transformed), vmax=np.amax(dG_transformed))
            # Changed units to microgals (multiply MKS unit by 10^8)
            self.norm = mcolor.Normalize(vmin=self.dmc['cbar_min'], vmax=self.dmc['cbar_max'])
        
        #self.m2.imshow(dG_transformed, cmap=cmap, norm=self.norm)
        # Changed units to microgals (multiply MKS unit by 10^8)
        if self.contours is None:
            self.m2.imshow(dG_transformed*float(pow(10,8)), cmap=cmap, norm=self.norm)
        else:
            map_x, map_y = self.m2(self.lons_1d, self.lats_1d)
            XX,YY = np.meshgrid(map_x, map_y)
            self.m2.contourf(XX, YY, dG_transformed*float(pow(10,8)), self.contours, cmap=cmap, norm=self.norm, extend='both')

        #-----------------------------------------------------------------------
        # Fig3 is the land/sea mask.
        #-----------------------------------------------------------------------
        """
        fig3 = mplt.figure(figsize=(mwi, mhi), dpi=plot_resolution)
        self.m3.ax = fig3.add_axes((0,0,1,1))
        #self.m3.fillcontinents(color='#000000', lake_color='#ffffff')
        dG_abs = np.fabs(dG_transformed)
        #print np.amin(dG_abs), np.amax(dG_abs), 1e6*np.amin(dG_abs), np.amax(dG_abs)*1e-1
        im = self.m3.imshow(dG_abs, cmap=mplt.get_cmap('gray_r'), norm=mcolor.Normalize(vmin=1e6*np.amin(dG_abs), vmax=np.amax(dG_abs)*1e-1, clip=True))
        
        #fig3.savefig('local/test_mask.png', format='png', dpi=plot_resolution)
        """
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
        """
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
        """
        
        # Clear all three figures
        fig1.clf()
        fig2.clf()
        #fig3.clf()
        mplt.close('all')
        gc.collect()
        
        #mask = Image.new("L", (w,h), 'black')
        # The final composited image.
        #return  Image.composite(im1, im2, mask)
        #return  Image.composite(im1, im2)
        return im2



