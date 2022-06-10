#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: ivandevert

GitHub repo: https://github.com/ivandevert/seisplot

channel naming convention: https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/
    letter order: 
        band code (sampling rate and response band of the instrument)
        instrument code (high/low gain, gravimeter, accelerometer, etc.)
        orientation code (Z/N/E, A/B/C, 1/2/3)
    obspy stream.select docs:
        https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.select.html

TO DO:
    add plot settings bar
    
    DPI is now hard coded to fix a bug with the plot not sizing correctly on 
    start (temp fix). A better solution should be implemented.
"""
import tkinter as tk
from tkinter import ttk, Grid
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use("TkAgg")

# efspy_module_parent_dir = "/Users/ivandevert/prog/efspy/resources/"

import sys
# sys.path.insert(1,efspy_module_parent_dir)
# sys.path.append("..")
# sys.path.append("../..")
import numpy as np
# import struct


import obspy
# from obspy import UTCDateTime, geodetics
from EFSpy_module import EFS
# from obspy.clients.fdsn import Client
import os
from pathlib import Path
# import csv
# from obspy import Catalog
# from obspy.core.event.event import Event
# from obspy.core.event import Origin
# from obspy.core.event import Magnitude
# import time
# import math
# from itertools import compress
# from datetime import date
# from datetime import datetime
# import timeit
import warnings
import copy
from tkinter import filedialog
import inspect
from configparser import ConfigParser


from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg #, NavigationToolbar2Tk
from matplotlib.figure import Figure

from scipy import signal
# from matplotlib.backend_bases import key_press_handler

# constants for program, don't change
LARGE_FONT= ("Verdana", 12)
DEF_GEOMETRY = [530, 1120]
DEF_CANVAS_GEOMETRY = [455, 845]

global CONFIG_FILENAME
CONFIG_FILENAME = 'config.ini'

global FIGNUM_SPECTROGRAM
global FIGNUM_PSD
global FIGNUM_RESPONSE
global FIGNUM_RECORD_SECTION

FIGNUM_SPECTROGRAM = 10
FIGNUM_PSD = 11
FIGNUM_RESPONSE = 12
FIGNUM_RECORD_SECTION = 13

def searchunsorted(A, v):
    # A is the array to search
    # v is the array to match into A
    
    sort = np.argsort(A)
    rank = np.searchsorted(A, v, sorter=sort)
    return sort[rank]

def get_geometry_string(geometry_list):
    """
    Returns formatted string for Tk geometry.
    Ex: [600, 800] returns '800x600'

    """
    lw = str(len(str(geometry_list[1])))
    lh = str(len(str(geometry_list[0])))
    return str('%'+lw+'ix%'+lh+'i') % (geometry_list[1], geometry_list[0])

def load_plugins(active_plugins):
    """
    Load plugin modules and store them in globals().
    
    This could be rewritten to be less work-aroundy

    """
    import importlib
    
    # store the heights of all plugins so window can be sized correctly
    global plugin_total_height
    plugin_total_height = 0
    
    # plugins path is just in the plugins/ subdirectory, add to system path
    pth = slash(str(Path().resolve())) + 'plugins/'
    sys.path.append(pth)
    
    plugins = []
    
    for el in os.listdir(pth):
        if Path(pth + el).is_dir() and el[0] != '_':
            if el in active_plugins:
                try:
                    globals()[el] = importlib.import_module(el+'.'+el, package=el)
                    plugins.append(el)
                    plugin_total_height += sys.modules['.'.join([el, el])].frame_height
                    print('Plugin ' + el + ' loaded successfully')
                except:
                    print("Error loading plugin " + el + ". Skipping.")
            else:
                # don't print status if plugin is sample_plugin
                if el not in ['sample_plugin']: 
                    print("Plugin " + el + " not active.")
    return plugins

def str2bool(string):
    """
    Turns a string 'True' or 'False' into a boolean True or False. Case insensitive.

    """
    S = string.lower()
    if S=='true':
        return True
    elif S=='false':
        return False
    else:
        raise ValueError("Acceptable string values are True or False (case insensitive)")

def init_config_file():
    """
    Generates a configuration file. This will only be run when config.ini does 
    not exist in the same directory.

    """
    config = ConfigParser()
    config.add_section('settings')
    config.set('settings', 'config_default_filepath', '~/')
    config.set('settings', 'config_p_wave_color', 'r')
    config.set('settings', 'config_s_wave_color', 'c')
    config.set('settings', 'config_debug_mode', 'False')
    config.set('settings', 'config_active_plugins', '')
    
    config.write(open('config.ini', 'w'))
    return

def load_config_handler():
    """
    Loads configuration settings from the configuration file. This is only used
    once per script execution.

    """
    # try: 
    #     debug_print('main')
    # except:
    #     print('DEBUG: main load_config (manual)')
    
    global CONFIG_FILENAME
    global config_default_filepath
    global config_p_wave_color
    global config_s_wave_color
    global config_debug_mode
    global active_plugins
    
    config = ConfigParser()
    
    if not os.path.isfile(CONFIG_FILENAME): 
        print("No config.ini file detected. Creating a new one.")
        init_config_file()

    try:
        config.read(CONFIG_FILENAME)
        load_config()
    except:
        print("Loading config.ini failed. Creating a new one and loading...", end='')
        try:
            init_config_file()
            load_config()
            print("Success")
        except:
            raise ValueError("Failed. You shouldn't see this error. Please contact me.")
    # print('DEBUGGING: ', config_debug_mode)
    return

def load_config():

    global CONFIG_FILENAME
    global config_default_filepath
    global config_p_wave_color
    global config_s_wave_color
    global config_debug_mode
    global active_plugins
    
    config = ConfigParser()

    config.read(CONFIG_FILENAME)
        
    config_default_filepath = config.get('settings', 'config_default_filepath').strip()
    config_p_wave_color = config.get('settings', 'config_p_wave_color').strip()
    config_s_wave_color = config.get('settings', 'config_s_wave_color').strip()
    config_debug_mode = str2bool(config.get('settings', 'config_debug_mode').strip())
    config_active_plugins = config.get('settings', 'config_active_plugins').split(',')
    active_plugins = [el.strip() for el in config_active_plugins]
    return

def save_config(setting_name, setting_value):
    """
    Saves a single configuration setting to the configuration file.

    Parameters
    ----------
    setting_name : string
        Name of the setting.
    setting_value : any (converts to string)
        Value to store.

    Returns
    -------
    None.

    """
    debug_print('main')
    global CONFIG_FILENAME
    
    config = ConfigParser()
    config.read(CONFIG_FILENAME)
    config.set('settings', setting_name, setting_value)
    # print("CHANGING: ", setting_name, " to ", str(setting_value))
    config.write(open(CONFIG_FILENAME, 'w'))
    return

def debug_print(parent_name):
    """
    If config_debug_mode is True, this will print the parent_name followed by
    the calling function's name. This is currently at the beginning of every
    function to allow for easy debugging.

    Parameters
    ----------
    parent_name : string
        Name of the parent class or function.

    Returns
    -------
    None.

    """
    global config_debug_mode
    global frame_tpf

    if config_debug_mode==True:
        print('DEBUG: ', parent_name, inspect.currentframe().f_back.f_code.co_name, end='')
        try:
            # print('\t fig: ', frame_tpf.fig.get_gid(), end='')
            print('\t nclick: ', frame_tpf.nclick, end='')
        except:
            print('', end='')

        print('\n', end='')

def slash(pth):
    """
    Appends a slash to the path string if necessary. Returns a slash if the 
    string is empty.

    Parameters
    ----------
    pth : string
        Path as a string.

    Returns
    -------
    string
        String where the last character is a slash.

    """
    if len(str(pth))>0:
        if str(pth)[-1]!='/':
            return str(pth)+'/'
        else:
            return str(pth)
    else:
        return '/'
    

def get_picks(tr):
    """
    Get all pick data from an input trace.
    
    POSSIBLE ISSUE: might not work with other file formats.
    TO DO: clean up and optimize if possible
    
    Parameters
    ----------
    tr : obspy.trace or efs.waveforms[n] (dict)
        Trace to get picks from.

    Returns
    -------
    pick_names, pick_times (relative to first sample), pick qualities

    """
    debug_print('main')
    # efs standard pick labels
    picklabels = ['pick1','pick2','pick3','pick4']
    # # pick time
    picks = np.array((),dtype=float)
    pick_names = np.array((),dtype=str)
    pick_qualities = np.array((),dtype=float)
        
    if type(tr)==obspy.core.trace.Trace:
        for pt in picklabels:
            if hasattr(tr.stats.pick_data,pt):
                picks = np.append(picks,tr.stats.pick_data[pt])
                pick_names = np.append(pick_names,tr.stats.pick_data[pt+"name"].strip())
                pick_qualities = np.append(pick_qualities,float(tr.stats.pick_data[pt+"q"]))
    
    # # for efs objects
    # elif type(tr)==dict:
    #     for pt in picklabels:
    #         if tr[pt+"name"] != '    ': # maybe use .strip() != '' or something
    #             picks = np.append(picks,tr[pt])
    #             pick_names = np.append(pick_names,tr[pt+"name"].strip())
    #             pick_qualities = np.append(pick_qualities,float(tr[pt+"q"]))
    return pick_names,picks,pick_qualities

def best_picks(pickn,pickt,pickq):
    """
    Returns one of each pick type corresponding to that pick type's highest
    quality pick.
    
    TO DO: pickq of 0.1 comes from my pick inserting function. This should be
    coded differently.

    Parameters
    ----------
    pickn : numpy array of strings
        Names of the picks (contents of each pick#name field).
    pickt : numpy array of floats
        Pick times.
    pickq : numpy array of floats
        Qualities of the picks.

    Returns
    -------
    pickn : numpy array of strings
        Pick names.
    pickt : numpy array of floats
        Pick times.
    pickq : numpy array of floats
        Pick qualities.

    """
    
    debug_print('main')
    # returns one of each of the pick types with the highest pick quality
    if len(pickn)>2:
        if len(pickn[pickn=='P'])>1:
            im = int(np.where(np.logical_and(pickn=='P',pickq==0.1))[0])
            pickn = np.delete(pickn,im)
            pickt = np.delete(pickt,im)
            pickq = np.delete(pickq,im)
        elif len(pickn[pickn=='S'])>1:
            im = int(np.where(np.logical_and(pickn=='S',pickq==0.1))[0])
            pickn = np.delete(pickn,im)
            pickt = np.delete(pickt,im)
            pickq = np.delete(pickq,im)
    
    return pickn,pickt,pickq

def get_files(path,filetype):
    """
    Return the full path of all files of a given filetype in a given directory.

    Parameters
    ----------
    path : string
        Directory to search for files in.
    filetype : string
        File extension to search for. Ex: '.efs', '.sac', etc.

    Returns
    -------
    list
        List of full paths of all files matching given filetype.

    """
    
    debug_print('main')
    filetype = filetype.lower()
    return [str(path)+el for el in os.listdir(path) if el[-len(filetype):].lower()==filetype]

def popup(title,message_str):
    """
    Generates a popup window with a given message. This does not format the
    string.

    Parameters
    ----------
    title : string
        Title of the popup window.
    message_str : string
        Message to display in the popup window.

    Returns
    -------
    None.

    """
    
    debug_print('main')
    win = tk.Toplevel()
    win.wm_title(title)

    l = tk.Label(win, text=message_str,anchor='w',justify='left',padx=20,pady=20)
    l.grid(row=0, column=0)

    b = ttk.Button(win, text="Okay", command=win.destroy)
    b.grid(row=1, column=0)
    
def popup_response(frame):
    """
    Plot the response of the filter. This should maybe be a function in TPF
    Does zerophase have an effect on the response? - don't think so
    """
    
    debug_print('main')
    tp = frame.pref_filter_type
    # freq = frame.pref_filter_freq
    f1 = frame.pref_filter_f1
    f2 = frame.pref_filter_f2
    f3 = frame.pref_filter_f3
    zerophase = frame.pref_filter_zerophase
    f_nyquist = frame.f_nyquist
    
    # print('DEBUG: ',tp,freq,f1,f2,f3,f_nyquist)
    
    if tp == 'highpass' or tp == 'lowpass':
        btype = tp
        F = f1/f_nyquist
    elif tp == 'bandpass':
        btype = 'band'
        F = [f1/f_nyquist, f2/f_nyquist]
        if F[1]>1: F[1]=1.0
    else:
        print('error showing response')
        return

    b, a = signal.iirfilter(f3,F,btype=btype,ftype='butter')
    w, h = signal.freqz(b,a,fs=f_nyquist*2,worN=512)
        
    fig = plt.figure(frame.fignum_response)
    plt.clf()
    ax = fig.add_subplot(111)
    ax.plot(w[1:-1], 20 * np.log10(np.abs(h[1:-1])))
    ax.set_title('Butterworth ' + tp + ' frequency response')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Amplitude (dB)')
    # ax.axis((0.01, f_nyquist*10, -100, 10))
    ax.grid(which='both', axis='both')
    plt.show()

class SeisPlot(tk.Tk):
    """
    Main frame that initializes and controls the other frames.
    
    """

    def __init__(self, *args, **kwargs):
        """
        Initialization function
        
        """

        global filepath
        global n_current_file
        
        global FIGNUM_SPECTROGRAM
        global FIGNUM_PSD
        global FIGNUM_RESPONSE
        global FIGNUM_RECORD_SECTION

        debug_print('SeisPlot')
        
                
        tk.Tk.__init__(self, *args, **kwargs)

        tk.Tk.wm_title(self, "seisplot client")
        
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand = True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        

        self.frames = {}

        for F in (StartPage, TracePlotFrame):

            frame = F(container, self)
            
            # print('frame '+str(F))
            self.frames[F] = frame

            frame.grid(row=0, column=0, sticky="nsew")
            
            
        self.show_frame(StartPage)

    def show_frame(self, cont):
        """
        Put a given frame on top of the gui

        Parameters
        ----------
        cont : tk.Frame
            Frame to show.

        Returns
        -------
        None.

        """
        debug_print('SeisPlot')

        show = True
        if cont==TracePlotFrame:
            if n_current_file >= 0:
                show = True
            else:
                show = False
        if show:
            frame = self.frames[cont]
            if cont==TracePlotFrame: frame.on_raise()
            frame.tkraise()
    
    def quit_program(self):
        """
        Holds code to run when gui is exited.
        

        """
        debug_print('SeisPlot')
        plt.close(FIGNUM_SPECTROGRAM)
        plt.close(FIGNUM_PSD)
        plt.close(FIGNUM_RESPONSE)
        plt.close(FIGNUM_RECORD_SECTION)
        self.destroy()
        self.quit()

class StartPage(tk.Frame):
    """
    Frame showing directory selection, plot, and quit buttons    
    
    """

    def __init__(self, parent, controller):
        """
        Initialize the StartPage frame.

        """
        debug_print('StartPage')
                
        global n_current_file
        global config_default_filepath
        global plugins
        
        n_current_file = -1
        self.controller = controller
        self.parent = parent
        self.plugins = plugins
        self.efs_files = []
        
        self.efs_parent_dir = config_default_filepath
        
        tk.Frame.__init__(self,parent)
        
        self.main_frame = tk.Frame(self)
        self.main_frame.pack()
        dir_select_label = tk.Label(self.main_frame, text="Directory selection", font=LARGE_FONT)
        

        self.plot_button = ttk.Button(self.main_frame, text="Plot traces",command=lambda: self.controller.show_frame(TracePlotFrame) )
        
        quit_button = ttk.Button(self.main_frame, text='Quit', command=lambda: self.controller.quit_program())
        
        label_efs_parent_dir = tk.Label(self.main_frame, text="EFS parent directory: ", font=LARGE_FONT)
        self.disp_efs_parent_dir = tk.Label(self.main_frame,text=self.efs_parent_dir)
        self.file_count_label = tk.Label(self.main_frame, text='null')
        button_change_efs_parent_dir = ttk.Button(self.main_frame,text='Change path',command=lambda: self.change_efs_parent_dir())
        
        
        dir_select_label.grid(row=0,columnspan=2)
        tk.Label(self.main_frame,text='').grid(row=2)
        label_efs_parent_dir.grid(row=3,columnspan=2)
        self.disp_efs_parent_dir.grid(row=4,columnspan=2)
        self.file_count_label.grid(row=5, columnspan=2)
        self.plot_button.grid(row=6,column=0,sticky='e')
        button_change_efs_parent_dir.grid(row=6,column=1,sticky='w')
        tk.Label(self.main_frame,text='').grid(row=7)
        quit_button.grid(row=8,columnspan=2)
        
        self.check_dir()
    
    def check_dir(self):
        """
        Checks if efs_parent_dir exists contains files and handles cases 
        accordingly.

        Raises
        ------
        ValueError
            Hopefully this won't be raised. If it is, idk what happened.

        TO DO: get_files is called here and in update_files() if nfiles > 0. 
            This could be slow with directories containing many files. Should 
            update this logic.

        """
        debug_print('StartPage')
        
        # case for the directory not existing
        if not os.path.isdir(self.efs_parent_dir): 
            self.file_count_label['text'] = 'Directory does not exist. Please choose a new directory.'
            self.file_count_label['fg'] = 'black'
            self.plot_button['state'] = 'disabled'
        
        # cases where the directory exists
        elif os.path.isdir(self.efs_parent_dir):
            
            nfiles = len(get_files(self.efs_parent_dir, '.efs'))
            
            # case where there are files in the directory
            if nfiles > 0:
                self.file_count_label['text'] = str(nfiles) + ' files found.'
                self.file_count_label['fg'] = 'green'
                self.plot_button['state'] = 'enabled'
                # print(self.efs_parent_dir)
                save_config('config_default_filepath', self.efs_parent_dir)
                self.update_files()
            
            # case where directory contains no files
            else:
                self.file_count_label['text'] = 'No files found in selected directory.'
                self.file_count_label['fg'] = 'red'
                self.plot_button['state'] = 'disabled'
        else:
            raise ValueError('Unknown error')
    
    def change_efs_parent_dir(self):
        """
        Generates a window to change the seismogram file parent directory and 
        calls check_dir() and populate_file_listbox() if it is not '/'
        
        Maybe populate_file_listbox() should be called in update_files()?

        Returns
        -------
        None.

        """
        debug_print('StartPage')
        folder_path = filedialog.askdirectory(initialdir=self.efs_parent_dir,title='Select the folder containing EFS files')
        
        if slash(folder_path) != '/':
            self.efs_parent_dir = slash(folder_path)
            self.check_dir()
            self.controller.frames[TracePlotFrame].populate_file_listbox()
        
    def update_files(self):
        """
        Updates the GUI objects relating to parent path, 

        Returns
        -------
        None.

        """
        debug_print('StartPage')
        global filepath
        global n_current_file
        
        # update project path
        self.disp_efs_parent_dir['text'] = self.efs_parent_dir
        
        # search dir for files and sort them     
        self.efs_files = get_files(self.efs_parent_dir,'.efs')
        self.efs_files.sort(reverse=False)
        
        # if directory has files, set globals n_current_file and filepath, 
        # otherwise set them to flag/null values
        if len(self.efs_files)>0:
            n_current_file = 0
            filepath = self.efs_files[n_current_file]
        else:
            n_current_file = -1
            filepath = 'null'
        

class TracePlotFrame(tk.Frame):
    """
    This is where everything is plotted
    
    """
    
    def __init__(self,parent,controller):
        """
        Initialize the TracePlotFrame frame
        
        The GUI is designed to have a single "parent parent" frame (self.main_frame), populated 
        with several gridded parent frames (labelled below), each with their own children and 
        layout. An exception is the canvas_frame, which is a tk.Canvas object and
        doesn't behave exactly the same.
        
        See https://tkdocs.com/tutorial/grid.html for info on grid and pack functions.        
        
        Also see the GUI layout figure in the /docs/ folder
        
        """
        
        debug_print('TracePlotFrame')
        global n_current_file
        global config_p_wave_color
        global config_s_wave_color
        
        global FIGNUM_SPECTROGRAM
        global FIGNUM_PSD
        global FIGNUM_RESPONSE
        global FIGNUM_RECORD_SECTION
        
        # store controller/parent in object so they can be used outside init
        self.controller = controller
        self.parent = parent
        
        # get the variables from parent frame
        self.efs_files = self.controller.frames[StartPage].efs_files
        self.efs_parent_dir = self.controller.frames[StartPage].efs_parent_dir
        self.plugins = self.controller.frames[StartPage].plugins
                
        self.NCOL = 14
                
        ### set default preferences. These should eventually be stored in the 
        # config.ini file
        
        # Filtering #
        self.pref_detrend = 1                # (0) none, (1) demean, (2) simple detrend
        self.pref_detrend_type = 'Demean'
        self.pref_detrend_order = 3
        self.pref_detrend_dspline = 1000
        self.pref_detrend_options = {}
        self.pref_filter_type = 'highpass'  # 'none', 'lowpass', 'highpass', 'bandpass'
        # self.pref_filter_f1 = 5           # default freq for lowpass and highpass
        self.pref_filter_f1 = 1        # default fmin for bandpass/highpass/lowpass
        self.pref_filter_f2 = 10       # default fmax for bandpass
        self.pref_filter_f3 = 4     # default corners for bandpass
        self.pref_filter_zerophase = False
        self.pref_id_filter_str = '*,*,*,*'
        self.pref_dist1_filter = 0
        self.pref_dist2_filter = 99999
        self.pref_trace_snr = 0
        
        self.pref_plot_picks = 'Real'       # T/F plot picks
        
        # convert pref into str. Can probably do this later and have GUI set to 'null' at first.
        # self.pref_detrend_str = self.pref_detrend_type
                
        # variables that hold values used for triggering
        self.manual_t_range = False
        self.changing_t_range = False
        self.nclick = 0
        
        # load other preferences
        self.p_c = config_p_wave_color
        self.s_c = config_s_wave_color
        
        ### program constants that really don't need to be changed, but I guess they could be
        self.scroll_percent = 10
        self.title_size = 8
        self.axes_label_size = 8
        self.axes_ticklabel_size = 6
        self.infobox_text_size = 6
        
        # set default figure numbers
        self.fignum_spectrogram = FIGNUM_SPECTROGRAM
        self.fignum_psd = FIGNUM_PSD
        self.fignum_response = FIGNUM_RESPONSE
        self.fignum_record_section = FIGNUM_RECORD_SECTION
        
        # GUI layout variables
        box_width = 6
        col_lb1_width = 12
        
        # initialize something? idk what this is for actually but the program breaks without it
        tk.Frame.__init__(self,parent)
        
#%% GUI Layout
        
        # "parent parent" frame
        self.main_frame = tk.Frame(self)
        self.main_frame.grid(row=0,column=0, sticky='nsew')
        
        #%% parent frames
        navbar_frame = tk.Frame(self.main_frame,padx=5,pady=5,borderwidth=5,relief='groove')
        file_nav_frame = tk.Frame(self.main_frame,padx=5,pady=5)
        # self.canvas_frame = tk.Frame(self.main_frame,padx=5,pady=5,borderwidth=2,relief='groove')
        self.canvas_frame = tk.Canvas(self.main_frame,borderwidth=2,relief='groove')
        trace_nav_frame = tk.Frame(self.main_frame,padx=5,pady=5)
        canvas_nav_frame = tk.Frame(self.main_frame,padx=5,pady=5)
        pref_frame = tk.Frame(self.main_frame,padx=5,pady=5,borderwidth=1,relief='groove')
        
        # grid parent frames
        navbar_frame.grid(column=0,row=0,columnspan=3,sticky='ew')
        file_nav_frame.grid(column=0,row=1,rowspan=3,sticky='nws')
        self.canvas_frame.grid(column=1,row=1,sticky='nsew')
        trace_nav_frame.grid(column=2,row=1,rowspan=3,sticky='nes')
        canvas_nav_frame.grid(column=1,row=2,sticky='s')
        pref_frame.grid(column=1,row=3,sticky='s')
        self.next_frame_row = 4
        self.plugin_colspan = 3
        
        #%% navigation frame
        # frame_label = tk.Label(navbar_frame, text="Trace plot",font=LARGE_FONT)
        # frame_label.pack(side='left')
        back_button = ttk.Button(navbar_frame, text='Back',command=lambda: controller.show_frame(StartPage))
        back_button.pack(side='left')
        
        save_button = ttk.Button(navbar_frame,text='Save figure',command=lambda: self.save_figure())
        save_button.pack(side='right')
        
        spectrogram_button = ttk.Button(navbar_frame,text='Plot spectrogram',command=lambda: self.plot_spectrogram())
        ppsd_button = ttk.Button(navbar_frame,text='Plot PSD',command=lambda: self.plot_ppsd())
        self.show_response_button = ttk.Button(navbar_frame,text='Filter response',command=lambda:popup_response(self))
        record_section_button = ttk.Button(navbar_frame, text='Plot record section', command=lambda: self.plot_record_section())
        
        tk.Label(navbar_frame,text=' ',width=10).pack(side='right')
        spectrogram_button.pack(side='right')
        ppsd_button.pack(side='right')
        record_section_button.pack(side='left')
        self.show_response_button.pack(side='right')
        ### end navigation frame
        
        #%% file nav frame
        file_nav_buttons_frame = tk.Frame(file_nav_frame)
        file_nav_buttons_frame.grid(row=0)
        button_previous_file = ttk.Button(file_nav_buttons_frame,text='Previous file',command=lambda: self.decrement_file())
        button_previous_file.pack(side='left')
        button_next_file = ttk.Button(file_nav_buttons_frame,text='Next file',command=lambda: self.increment_file())
        button_next_file.pack(side='right')
        
        choose_file_frame = tk.Frame(file_nav_frame)
        choose_file_frame.grid(row=1,sticky='ns')
        choose_file_label = tk.Label(choose_file_frame,text='Select a file: ')
        choose_file_label.pack(side='top')
        
        file_listbox_frame = tk.Frame(file_nav_frame)
        file_listbox_frame.grid(row=2,sticky='ns')
        self.file_listbox = tk.Listbox(file_listbox_frame,height=23,exportselection=False)
        self.file_listbox.pack(side='left',fill='both')
        file_scrollbar = tk.Scrollbar(file_listbox_frame)
        file_scrollbar.pack(side='right',fill='both')
        
        self.file_listbox.config(yscrollcommand=file_scrollbar.set)
        file_scrollbar.config(command=self.file_listbox.yview)
        
        self.file_listbox_label = tk.Label(file_nav_frame,text='file counter')
        self.file_listbox_label.grid(row=3,sticky='n')
        ### end file nav frame
        
        
        
        
        #%% canvas frame
        # old
        self.fig = Figure( figsize=np.flip(DEF_CANVAS_GEOMETRY)/75, dpi=75, frameon=False, layout='tight')
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.grid(row=0, column=0, sticky='nsew')
        self.canvas_frame.rowconfigure(0, weight=10)
        self.canvas_frame.columnconfigure(0, weight=10)
        # self.canvas_frame.config(width=DEF_CANVAS_GEOMETRY[1], height=DEF_CANVAS_GEOMETRY[0])
        # self.canvas_frame.itemconfigure(self.canvas, width=DEF_CANVAS_GEOMETRY[1], height=DEF_CANVAS_GEOMETRY[0])
        ### end canvas frame
        
        #%% canvas navigation frame
        self.zoom_toggle_button = ttk.Button(canvas_nav_frame,text='Change time range', command=lambda: self.on_zoom_toggle_click())
        self.zoom_toggle_button.pack(side='left')
        self.zoom_status_label = tk.Label(canvas_nav_frame,text='',fg='red',width=20)
        self.zoom_status_label.pack(side='left')
        
        # scroll frame
        self.scroll_frame = tk.Frame(canvas_nav_frame)
        
        self.scroll_label = tk.Label(self.scroll_frame,text='Scroll')
        self.scroll_left_button = ttk.Button(self.scroll_frame,text='<-',command=lambda: self.scroll(-1),width=2)
        self.scroll_right_button = ttk.Button(self.scroll_frame,text='->',command=lambda: self.scroll(1),width=2)
        self.scroll_left_button.pack(side='left')
        self.scroll_label.pack(side='left')
        self.scroll_right_button.pack(side='left')
        # end scroll frame
        
        ### plot picks
        # plot_picks_values = ['None','Real']
        plot_picks_values = ["None", "Real"]
        # plot_picks_label = tk.Label(canvas_nav_frame,text='Plot picks: ',width=col_lb1_width,anchor='w')
        # plot_picks_label.grid(column=0,row=2,sticky='w')
        current_pref_plot_picks = tk.StringVar()
        self.plot_picks_box = ttk.Combobox(canvas_nav_frame,textvariable=current_pref_plot_picks,values=plot_picks_values,state='readonly',width=10)
        self.plot_picks_box.pack(side='right')
        self.plot_picks_box.set(self.pref_plot_picks)
        ### end canvas navigation frame
        
        self.scroll_frame.pack(side='right')
        
        #%% trace navigation frame
        trace_nav_buttons_frame = tk.Frame(trace_nav_frame)
        trace_nav_buttons_frame.grid(row=0)
        next_trace_button = ttk.Button(trace_nav_buttons_frame,text="Next trace",command= lambda:self.increment_trace())
        next_trace_button.pack(side='right')
        previous_trace_button = ttk.Button(trace_nav_buttons_frame,text="Previous trace",command= lambda:self.decrement_trace())
        previous_trace_button.pack(side='left')
        
        choose_trace_frame = tk.Frame(trace_nav_frame)
        choose_trace_frame.grid(row=1,sticky='ns')
        choose_trace_label = tk.Label(choose_trace_frame,text='Select a trace: ')
        choose_trace_label.pack(side='top')
        
        trace_listbox_frame = tk.Frame(trace_nav_frame)
        trace_listbox_frame.grid(row=2,sticky='n')
        self.trace_listbox = tk.Listbox(trace_listbox_frame,height=17,exportselection=False)
        self.trace_listbox.pack(side='left',fill='both')
        trace_scrollbar = tk.Scrollbar(trace_listbox_frame)
        trace_scrollbar.pack(side='right',fill='both')
        self.trace_listbox.config(yscrollcommand=trace_scrollbar.set)
        trace_scrollbar.config(command=self.trace_listbox.yview)
        
        self.trace_listbox_label = tk.Label(trace_nav_frame,text='trace counter')
        self.trace_listbox_label.grid(row=3,sticky='n')
        
        # trace filtering
        trace_filter_frame = tk.Frame(trace_nav_frame)
        trace_filter_frame.grid(row=4,sticky='ns')
        
        trace_filter_frame_label = tk.Label(trace_filter_frame,text='Only show traces where: ',anchor='w')
        trace_filter_id_label = tk.Label(trace_filter_frame,text='ID: ',anchor='w')
        trace_filter_distance_label = tk.Label(trace_filter_frame,text='Distance (km): ',anchor='w')
        trace_filter_to_label = tk.Label(trace_filter_frame,text='to',width=1)
        trace_snr_filter_label = tk.Label(trace_filter_frame, text='Minimum SNR: ', anchor='w')
        
        component_var = tk.StringVar()
        dist1_var = tk.StringVar()
        dist2_var = tk.StringVar()
        snr_var = tk.StringVar()
        self.trace_filter_id_entry = ttk.Entry(trace_filter_frame,textvariable=component_var,width=8)
        self.trace_filter_dist1_entry = ttk.Entry(trace_filter_frame,textvariable=dist1_var,width=2)
        self.trace_filter_dist2_entry = ttk.Entry(trace_filter_frame,textvariable=dist2_var,width=5)
        self.trace_filter_id_help_button = ttk.Button(trace_filter_frame,text='?',width=0.5,command=lambda: self.filter_help_popup())
        self.trace_snr_filter_entry = ttk.Entry(trace_filter_frame, textvariable=snr_var, width=3)
        
        trace_filter_frame_label.grid(row=0,column=0,columnspan=4,sticky='w')
        trace_filter_id_label.grid(row=1,column=0,sticky='w')
        trace_filter_distance_label.grid(row=2,column=0,sticky='w')
        trace_filter_to_label.grid(row=2,column=2)
        trace_snr_filter_label.grid(row=3, column=0, sticky='w')
        
        self.trace_filter_id_entry.grid(row=1,column=1,columnspan=2,sticky='w')
        self.trace_filter_id_help_button.grid(row=1,column=3,sticky='e')
        self.trace_filter_dist1_entry.grid(row=2,column=1)
        self.trace_filter_dist2_entry.grid(row=2,column=3)
        self.trace_snr_filter_entry.grid(row=3, column=1)
        
        self.trace_filter_id_entry.insert(0,self.pref_id_filter_str)
        self.trace_filter_dist1_entry.insert(0,str(self.pref_dist1_filter))
        self.trace_filter_dist2_entry.insert(0,str(self.pref_dist2_filter))
        self.trace_snr_filter_entry.insert(0, str(self.pref_trace_snr))
        # end trace filtering
        
        
        # # sorting MAYBE COMING SOON idk
        # self.trace_order_text = 'NULL'
        # self.trace_sort_ascend = True
        # selected = tk.StringVar()
        # trace_sort_frame = tk.Frame(trace_nav_frame)
        # trace_sort_frame.grid(row=4)
        
        # self.trace_sort_order_button = ttk.Button(trace_sort_frame,text='Sorting by ('+self.trace_order_text+'): ',command=lambda: self.on_change_trace_order_press())
        # self.radio_sort_abc = ttk.Radiobutton(trace_sort_frame, text='Alphabetical', value='abc', variable=selected)
        # self.radio_sort_distance = ttk.Radiobutton(trace_sort_frame, text='Distance', value='distance', variable=selected)
        
        # self.trace_sort_order_button.grid(row=0)
        # # end sorting
        ### end trace navigation frame
        
        #%% pref frame
        
        # detrend row row=0
        detrend_values = ['None', 'Simple', 'Demean', 'Linear', 'Polynomial', 'Spline']
        current_pref_detrend = tk.StringVar()
        detrend_label = ttk.Label(pref_frame,text='Detrend option: ',width=col_lb1_width,anchor='w')
        self.detrend_box = ttk.Combobox(pref_frame, textvariable=current_pref_detrend,values=detrend_values,state='readonly',width=box_width)
        self.detrend_box.set(self.pref_detrend_type)
        self.detrend1_entry = ttk.Entry(pref_frame, textvariable=tk.StringVar(), width=4)
        self.detrend2_entry = ttk.Entry(pref_frame, textvariable=tk.StringVar(), width=4)
        self.detrend1_label = ttk.Label(pref_frame, text='Order: ', anchor='w')
        self.detrend2_label = ttk.Label(pref_frame, text='dspline: ', anchor='w')
        
        detrend_label.grid(column=0,row=0,sticky='w')
        self.detrend_box.grid(column=1,row=0,sticky='w')
        self.detrend1_label.grid(column=2, row=0, sticky='w')
        self.detrend1_entry.grid(column=3, row=0, sticky='w')
        self.detrend2_label.grid(column=4, row=0, sticky='w')
        self.detrend2_entry.grid(column=5, row=0, sticky='w')
        
        self.detrend1_entry.insert(0, self.pref_detrend_order)
        self.detrend2_entry.insert(0, self.pref_detrend_dspline)
        self.detrend1_entry['state'] = 'disabled'
        self.detrend2_entry['state'] = 'disabled'
        # end detrend row
        
        # filter row row=1
        filter_type_values = ['none','highpass','lowpass','bandpass']
        current_pref_filter_type = tk.StringVar()
        filter_type_label = tk.Label(pref_frame,text="Filter type: ",width=col_lb1_width,anchor='w')
        self.filter_type_box = ttk.Combobox(pref_frame,textvariable=current_pref_filter_type,values=filter_type_values,state='readonly',width=box_width)
        self.filter_type_box.set(self.pref_filter_type)
        
        self.f1_label = tk.Label(pref_frame,text='null',anchor='e')
        self.f2_label = tk.Label(pref_frame,text='null',anchor='e') # ,width=col_lb1_width
        self.f3_label = tk.Label(pref_frame,text='null',anchor='e')
        
        self.f1_entry = ttk.Spinbox(pref_frame,width=3,from_=0.001,to=100.0)
        self.f2_entry = ttk.Spinbox(pref_frame,width=3,from_=0.0,to=100.0)
        self.f3_entry = ttk.Spinbox(pref_frame,width=3,from_=0.0,to=100.0)
        self.zerophase_checkbutton_value = tk.IntVar()
        self.zerophase_checkbutton = ttk.Checkbutton(pref_frame, onvalue=1, offvalue=0, command=lambda: self.on_zerophase_checkbutton_toggle(self.zerophase_checkbutton_value), text=' Zerophase')
        
        self.f1_entry.insert(0,str(self.pref_filter_f1))
        self.f2_entry.insert(0,str(self.pref_filter_f2))
        self.f3_entry.insert(0,str(self.pref_filter_f3))
        
        # grid filter row
        filter_type_label.grid(column=0,row=1,sticky='w')
        self.filter_type_box.grid(column=1,row=1,sticky='w')
        self.f1_label.grid(column=2,row=1,sticky='w')
        self.f1_entry.grid(column=3,row=1,sticky='w')
        self.f2_label.grid(column=4,row=1,sticky='w')
        self.f2_entry.grid(column=5,row=1,sticky='w')
        self.f3_label.grid(column=6,row=1,sticky='w')
        self.f3_entry.grid(column=7,row=1,sticky='w')
        self.zerophase_checkbutton.grid(column=8, row=1, sticky='e')
        # self.show_response_button.grid(column=8,row=1,sticky='e')
        
        pref_frame.rowconfigure(1, weight=1)
        # end filter row
        
        self.on_filter_option_change(True)
        
        self.plugins_grid()
        
        #%% these lines allow canvas_frame to change shape on window resize
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)
        self.main_frame.rowconfigure(1, weight=1)
        self.main_frame.columnconfigure(1, weight=1)
        # self.canvas.draw()
        
        # end plot picks row
        #%% end pref frame
        
#%%
        
        self.n_current_trace = 0
        self.t0 = 0.0
        self.tf = -99999.
        
        # let the plugins initialize
        self.plugins_on_traceplotframe_init()
        
        # if n_current_file isn't a -1 flag, 
        if n_current_file >= 0:
            self.populate_file_listbox()
        
        # this used to be a repeat, since on_file_change() is called in populate_file_listbox()
        # print(n_current_file)
        # if n_current_file >= 0:
        #     self.on_file_change()
    
    def on_raise(self):
        self.refresh()
        return
#%% PLUGIN HANDLERS
    
    def plugins_grid(self):
        """
        Grid the plugins

        """
        debug_print('TracePlotFrame n='+str(len(plugins)))
        nplugins = len(plugins)
        
        self.plugin_frames = [[]] * nplugins
        
        for ii, pl in enumerate(self.plugins):
            
            try:
                # initialize a frame and grid it
                self.plugin_frames[ii] = tk.Frame(self.main_frame,padx=5,pady=5,borderwidth=5,relief='groove')
                # self.plugin_frames[ii] = tk.Frame(self.main_frame,padx=5,pady=5)
                self.plugin_frames[ii].grid(column=0, row=self.next_frame_row, columnspan=self.plugin_colspan, sticky='ew')
    
                # send frame to plugin's on_grid() function to fill in
                sys.modules['.'.join([pl, pl])].on_grid(self, self.plugin_frames[ii])
                
                self.next_frame_row += 1
            except:
                warnings.warn("Problem gridding plugin " + str(pl) + ". Skipping.")
            
        return
    
#%% FILE NAVIGATION FUNCTIONS
    
    def populate_file_listbox(self):
        """
        Populates the file listbox with filenames in the efs_files list.

        """
        debug_print('TracePlotFrame')
        global n_current_file
        
        n_current_file = 0
        self.efs_files = self.controller.frames[StartPage].efs_files
        
        # filenames = [el.split('/')[-1].split('.')[0] for el in self.efs_files] # this lists the filenames without the extension
        filenames = [el.split('/')[-1] for el in self.efs_files]
        self.nfiles = len(filenames)
        self.file_listbox.delete(0,'end')
        for el in filenames:
            self.file_listbox.insert('end',el)
        self.file_listbox.select_set(0)
        self.file_listbox_label['text'] = str(n_current_file+1)+" of "+str(len(self.efs_files))+" files"
        self.on_file_change()
    
    def increment_file(self):
        """
        Increment n_current_file and change accordingly

        """
        debug_print('TracePlotFrame')
        global n_current_file
        if n_current_file < len(self.efs_files)-1:
            n_current_file += 1
            self.file_listbox.selection_clear(0,'end')
            self.file_listbox.select_set(n_current_file)
            self.on_file_change()
    
    def decrement_file(self):
        """
        Decrement n_current_file and change accordingly

        """
        debug_print('TracePlotFrame')
        global n_current_file
        if n_current_file > 0:
            n_current_file -= 1
            self.file_listbox.selection_clear(0,'end')
            self.file_listbox.select_set(n_current_file)
            self.on_file_change()
    
    def on_file_change(self):
        debug_print('TracePlotFrame')
        global n_current_file
        
        self.new_file = True

        self.efs_files = self.controller.frames[StartPage].efs_files
        
        n_current_file = self.file_listbox.curselection()[0]
        self.n_current_trace = 0
        
        self.reset_zoom()
        
        self.load_file()
        self.fill_fields()
        self.filter_trace_list()
        
        self.file_listbox_label['text'] = str(n_current_file+1)+" of "+str(len(self.efs_files))+" files"
        
        self.refresh_st()
        self.populate_trace_listbox(self.st_pref_subset)
    
    def load_file(self):
        debug_print('TracePlotFrame')
        global n_current_file
        global efs
        
        # this will soon be an issue (why?)
        self.filepath = self.efs_files[n_current_file]
        # print('ncurr: ',n_current_file)
        self.efs = EFS(self.filepath,np.float32,np.int32)
        self.efs_pref = copy.deepcopy(self.efs)
        efs = self.efs
        self.st = self.efs_pref.to_obspy(keep_pkdata=True)
        self.st.sort()
        self.st_subset = self.st.copy()
        self.st_pref = self.st_subset.copy()
        self.st_pref_subset = self.st_pref.copy()
        # self.tr_show_bool = np.ones(len(self.st), dtype=bool)
        # self.tr_ids_order = [el.get_id() for el in self.st.traces]
        # if len(self.tr_ids_order) != len(set(self.tr_ids_order)): raise ValueError("Not all EFS ids unique")
        
        
#%% TRACE NAVIGATION FUNCTIONS
    
    def populate_trace_listbox(self, st):
        """
        Populates the trace listbox with traces in the st_subset Stream object.

        """
        debug_print('TracePlotFrame')
        trace_ids = [str(el).split('|')[0].strip() for el in st.traces]

        # change selection to the first one if the traces have changed
        if len(trace_ids) != self.ntr:
            self.n_current_trace = 0

        self.ntr = len(trace_ids)
        self.trace_listbox.delete(0,'end')
        
        if self.ntr > 0:
            for el in trace_ids:
                self.trace_listbox.insert('end',el)
        else:
            self.trace_listbox.insert('end','No traces available')
        self.trace_listbox.select_set(self.n_current_trace)
        # self.trace_listbox.select_set(self.n_current_trace)
        
    def increment_trace(self):
        debug_print('TracePlotFrame')
        if self.n_current_trace<self.ntr-1: 
            self.n_current_trace += 1
            self.trace_listbox.selection_clear(0,'end')
            self.trace_listbox.select_set(self.n_current_trace)
            self.on_trace_change()
            self.refresh()
        
    def decrement_trace(self):
        debug_print('TracePlotFrame')
        if self.n_current_trace>0: 
            self.n_current_trace -= 1
            self.trace_listbox.selection_clear(0,'end')
            self.trace_listbox.select_set(self.n_current_trace)
            self.on_trace_change()
            self.refresh()

    def on_trace_change(self):
        debug_print('TracePlotFrame')
        self.manual_t_range = False
        nsel = self.trace_listbox.curselection()[0]
        self.n_current_trace = nsel
        self.reset_zoom()
        self.refresh()
    
    # might use modified version of this in future (not working now)
    # def order_traces(self, st):
    #     tr_ids = [el.get_id() for el in st.traces]
    #     ids_ordered = searchunsorted(self.tr_ids_order, tr_ids)
    #     st_ordered = st.select(id=ids_ordered)
    #     return st_ordered
    
    def on_distance_or_id_filter_change(self):
        debug_print('TracePlotFrame')
        
        snr = float(self.trace_snr_filter_entry.get())
        # print("SNR----------->: ", snr)
        
        self.filter_trace_list()
        
        # refresh the trace listbox and stream, and set the current trace to 0
        # if snr option is set, skip populate_trace_listbox() (because it will be called later)
        if snr < 0.0001:
            self.populate_trace_listbox(self.st_pref_subset)
        # self.n_current_trace = 0
        self.on_trace_change()
        self.refresh_st()
        
        return
    
    def filter_trace_list(self):
        """
        Filter traces out of the trace list. Right now, you can filter by distance
        or trace id contents.

        """
        debug_print('TracePlotFrame')
        
        # get the entered user values
        id_str = self.trace_filter_id_entry.get()
        distance1 = float(self.trace_filter_dist1_entry.get())
        distance2 = float(self.trace_filter_dist2_entry.get())
        # snr = float(self.trace_snr_filter_entry.get())
        
        # self.tr_show_bool
        
        # if nothing changed, return.
        if id_str == self.pref_id_filter_str and distance1 == self.pref_dist1_filter and distance2 == self.pref_dist2_filter and self.new_file==False: 
            # print('unchanged')
            return
        
        # otherwise, select only the traces fitting the parameters
        else:
            try:
                self.pref_id_filter_str = id_str
                self.pref_dist1_filter = distance1
                self.pref_dist2_filter = distance2
                
                network_str,station_str,location_str,channel_str = id_str.lower().split(',')
                self.st_subset = self.st.select(network=network_str,station=station_str,location=location_str,channel=channel_str)
                
                for tr in self.st_subset:
                    deldist = tr.stats.station_data['deldist']
                    if deldist < distance1 or deldist > distance2:
                        self.st_subset.remove(tr)                        
                    
            except Exception as err:
                print("Invalid string entry. Please try again (click ? for help)")
                print(err)
                return
        
        
        # sort traces back in order
        self.st_subset.sort()   
        self.ntr = len(self.st_subset)
        self.n_current_trace = 0
        return
    
    def on_snr_filter_change(self):
        debug_print('TracePlotFrame')
        
        self.filter_trace_list_snr()
        # self.n_current_trace = 0
        self.populate_trace_listbox(self.st_pref_subset)
        self.refresh_st()
        return
    
    def on_detrend_option_change(self):
        debug_print('TracePlotFrame')
        pref_detrend_type = self.detrend_box.get()
        pref_detrend_order = self.detrend1_entry.get()
        pref_detrend_dspline = self.detrend2_entry.get()
        
        # only update values if anything has changed
        if np.logical_or(pref_detrend_type != self.pref_detrend_type, np.logical_or(pref_detrend_order != self.pref_detrend_order, pref_detrend_dspline != self.pref_detrend_dspline)):
            self.pref_detrend_type = pref_detrend_type
            self.pref_detrend_order = pref_detrend_order
            self.pref_detrend_dspline = pref_detrend_dspline
            
            if self.pref_detrend_type == 'Polynomial':
                self.pref_detrend_options = {'order':int(self.pref_detrend_order)}
                self.detrend1_entry['state'] = 'enabled'
                self.detrend2_entry['state'] = 'disabled'
                self.detrend1_label['text'] = 'Order'
                self.detrend2_label['text'] = '----'
    
            elif self.pref_detrend_type == 'Spline':
                self.pref_detrend_options = {'order':int(self.pref_detrend_order), 'dspline':int(self.pref_detrend_dspline)}
                self.detrend1_entry['state'] = 'enabled'
                self.detrend2_entry['state'] = 'enabled'
                self.detrend1_label['text'] = 'Order'
                self.detrend2_label['text'] = 'dspline'
            else:
                self.pref_detrend_options = {}
                self.detrend1_entry['state'] = 'disabled'
                self.detrend2_entry['state'] = 'disabled'
                self.detrend1_label['text'] = '----'
                self.detrend2_label['text'] = '----'
    
            
            self.refresh_st()
        return
    
    def on_filter_option_change(self, init):
        debug_print('TracePlotFrame')
        
        pref_filter_type = self.filter_type_box.get()
        f1_value = float(self.f1_entry.get())
        f2_value = float(self.f2_entry.get())
        f3_value = int(self.f3_entry.get())
        pref_filter_zerophase = bool(self.zerophase_checkbutton_value.get())
        
        # if init, self. values don't exist yet, so cont is an intermediate control variable
        if init==True:
            cont = True
        # check if settings have changed
        elif int(pref_filter_type != self.pref_filter_type) + int(f1_value != self.f1_value) + int(f2_value != self.f2_value) + int(f3_value != self.f3_value) + int(pref_filter_zerophase != self.pref_filter_zerophase) > 0:
            cont = True
        else:
            cont = False
        
        if cont == True:
            self.pref_filter_type = pref_filter_type
            self.f1_value = f1_value
            self.f2_value = f2_value
            self.f3_value = f3_value
            self.pref_filter_zerophase = pref_filter_zerophase
            
            # get values, update labels and entry fields
            if self.pref_filter_type=='none':
                self.f1_label['text'] = '----'
                self.f2_label['text'] = '----'
                self.f3_label['text'] = '----'
                self.f1_entry['state'] = 'disabled'
                self.f2_entry['state'] = 'disabled'
                self.f3_entry['state'] = 'disabled'
                self.zerophase_checkbutton['state'] = 'disabled'
                self.show_response_button['state'] = 'disabled'
            elif self.pref_filter_type=='highpass' or self.pref_filter_type=='lowpass':
                self.pref_filter_f1 = self.f1_value
                
                self.f2_entry['state'] = 'disabled'
                self.f3_entry['state'] = 'disabled'
                self.f1_label['text'] = 'Frequency: '
                self.f2_label['text'] = '----'
                self.f3_label['text'] = '----'
                self.f1_entry['state'] = 'enabled'
                self.f2_entry['state'] = 'disabled'
                self.f3_entry['state'] = 'disabled'
                self.zerophase_checkbutton['state'] = 'enabled'
                self.show_response_button['state'] = 'enabled'
    
            elif self.pref_filter_type=='bandpass':
                self.pref_filter_f1 = self.f1_value
                self.pref_filter_f2 = self.f2_value
                self.pref_filter_f3 = self.f3_value
                
                self.f1_label['text'] = 'F min: '
                self.f2_label['text'] = 'F max: '
                self.f3_label['text'] = 'Corners: '
                self.f1_entry['state'] = 'enabled'
                self.f2_entry['state'] = 'enabled'
                self.f3_entry['state'] = 'enabled'
                self.zerophase_checkbutton['state'] = 'enabled'
                self.show_response_button['state'] = 'enabled'
                
            # if this isn't being run on the frame's __init__ function, call refresh_st()
            if init != True: self.refresh_st()
        
        return
    
    def filter_trace_list_snr(self):
        """
        Filter out traces based on SNR preference. Calculates on st_pref stream. Unusable in current state
        Issue: SNR should be calculated after filtering. filter_trace_list happens before filtering.
        idk what the best route is here. Maybe filter stream upon first load, calculate snr, and just use
        those values?

        """
        debug_print('TracePlotFrame')
        
        snr = float(self.trace_snr_filter_entry.get())
        self.st_pref_subset = self.st_pref.copy()
        
        # if snr preference hasn't changed, don't do anything
        if snr==0.0:
            # print('snr not set')
            return
        else:
            self.pref_trace_snr = snr
            self.n_current_trace = 0
            # print('snr set to ', snr)
            for tr in self.st_pref_subset:
                
                if hasattr(tr.stats.pick_data,'tdif'):
                    tdif = tr.stats.pick_data['tdif']
                else:
                    tdif = 0
                
                t = tr.times() + tdif
                
                t0 = t[0]
                tf = t[-1]
                
                delta_t = tf - t0
                time_buffer = 3*(self.perc_time_buffer/100)*delta_t
                # print('t window: ',self.t0+time_buffer,self.tf-time_buffer)
                dwin = tr.data[np.logical_and((t>=t0+time_buffer),(t<=tf-time_buffer))]
                
                trmax = np.max(np.abs(dwin))
                trmean = np.mean(np.abs(dwin))
                
                if trmax/trmean < snr:
                    self.st_pref_subset.remove(tr)
            self.ntr = len(self.st_pref_subset)
        return
            # self.populate_trace_listbox(self.st_pref_subset)
        
        
        
    # this might be used in the future so I'll leave it for now
    # def on_change_trace_order_press(self):
        
    #     if self.trace_sort_ascend:
    #         self.trace_sort_ascend = False
    #         self.trace_order_text = 'descend'
    #     else:
    #         self.trace_sort_ascend = True
    #         self.trace_order_text = 'ascend'
    #     self.trace_sort_order_button['text'] = 'Sorting by ('+self.trace_order_text+'): '
    
#%% POPUP WINDOW FUNCTIONS
    
    def plot_spectrogram(self):
        debug_print('TracePlotFrame')
        fig = plt.figure(self.fignum_spectrogram)
        plt.clf()
        ax = fig.add_subplot(111)
        tr = self.tr_pref.copy()
        tr.spectrogram(axes=ax, show=False)
        mappable = ax.images[0]
        ax.set_xlabel('Time (s)')
        ax.set_ylabel("Frequency (Hz)")
        ax.set_title(self.title_str)
        cb = plt.colorbar(mappable=mappable, ax=ax)
        cb.set_label("Intensity")
        plt.tight_layout()
        plt.show()
        return
    
    def plot_ppsd(self):
        debug_print('TracePlotFrame')
        fig = plt.figure(self.fignum_psd)
        plt.clf()
        # ax = fig.add_subplot(111)
        tr = self.tr_pref
        d = tr.data
        Pxx, freqs = plt.psd(d,Fs=tr.stats.sampling_rate, figure=fig)
        plt.show()
        return
    
    def plot_record_section(self):
        debug_print('TracePlotFrame')
        
        try:
            fig = plt.figure(self.fignum_record_section)
            plt.clf()
            if hasattr(self.st[0].stats,'distance'):
                st_record = self.st_pref_subset
            else:
                st_record = self.st_pref_subset
                for ii in np.arange(len(st_record)):
                    st_record[ii].stats.distance = st_record[ii].stats.station_data.deldist * 1000
            st_record.plot(type='section', fig=fig, orientation='horizontal')
            plt.show()
        except:
            print('error with record section')
        return
    
    def update_external_figures(self):
        debug_print('TracePlotFrame')
        fn = plt.get_fignums()
        if self.fignum_spectrogram in fn: 
            self.plot_spectrogram()
            # print('spectrogram updated')
        if self.fignum_psd in fn: 
            self.plot_ppsd()
            # print('psd updated')
        if self.fignum_response in fn: 
            popup_response(self)
            # print('response updated')
        if self.fignum_record_section in fn: 
            self.plot_record_section()
            # print('record section updated')
        return
        
    def filter_help_popup(self):
        """
        Displays a simple help popup for filtering the trace list

        """
        debug_print('TracePlotFrame')
        lines = ['Enter Unix-style, comma-separated strings to filter out certain traces. ',
                 '* is a wildcard for one or more characters',
                 '? is a wildcard for exactly one character',
                 'bracketed [] characters match on any character inside the brackets',
                 " ",
                 'Format: \"NET,STA,LOC,CHA\"',
                 'where each of the above comma-separated values represents a search string',
                 'for the corresponding property.',
                 ' ',
                 'Example: \"*,*,*,H?[ew]\"',
                 '\t- All networks, stations, and locations are kept',
                 '\t- Channels whose first character is H, second character is anything,',
                 '\t  and third character is either E or W are kept',
                 ' ',
                 'See documentation for obspy.core.stream.Stream.select for more info.']
        
        message_str = '\n'.join(lines)
        popup('Station ID filtering',message_str)
    
#%% PLOT MANIPULATION FUNCTIONS
    
    def on_zoom_toggle_click(self):
        debug_print('TracePlotFrame')
        
        # If t range is already being manually set, reset the zoom
        if self.manual_t_range:
            self.reset_zoom()
            
        # Otherwise, start zoom process
        else:
            self.manual_t_range = True
            self.changing_t_range = True
            self.nclick = 0
            self.zoom_toggle_button['text'] = 'Reset zoom'
            self.zoom_status_label['text'] = 'Select t1'
            
        self.refresh()
    
    def reset_zoom(self):
        debug_print('TracePlotFrame')
        self.tf = -99999.0
        self.t0 = 0.0
        self.zoom_status_label['text'] = ''
        self.manual_t_range = False
        self.changing_t_range = False
        self.zoom_toggle_button['text'] = 'Change time range'
    
    
    def on_canvas_click(self, event):
        debug_print('TracePlotFrame')
        # print(self.nclick,[event.xdata,event.ydata])
        # print(event)
        # if on manual zoom
        if self.manual_t_range:
            if event.inaxes is not None:
                self.nclick += 1
                # print("THIS, nclick = ", self.nclick)
                if self.nclick==1:
                    self.t1 = event.xdata
                    # self.nclick += 1
                    self.zoom_status_label['text'] = 'Select t2'
                elif self.nclick==2:
                    self.t2 = event.xdata
                    trange = [self.t1,self.t2]
                    trange.sort()
                    # print('trange: ',trange)
                    self.t0 = trange[0]
                    self.tf = trange[1]
                    self.nclick = 0
                    self.zoom_status_label['text'] = 'Select t1'
                    self.changing_t_range = False
                    self.refresh()
                
            else:
                print('clicked outside of axes bounds but inside plot window')
    
    def scroll(self,direction):
        """
        Scroll the plot window left or right

        """
        debug_print('TracePlotFrame')
        #        initial +  (-1/1)  * scroll magnitude       * window width
        t0_new = self.t0 + direction*(self.scroll_percent/100)*(self.tf-self.t0)
        tf_new = self.tf + direction*(self.scroll_percent/100)*(self.tf-self.t0)
        
        # if trying to scroll too far left
        if t0_new <= self.tmin:
            t_shift = direction*(self.t0-self.tmin)
        
        # if new t0, tf are between tmin, tmax
        elif t0_new > self.tmin and tf_new < self.tmax:
            t_shift = direction * (self.scroll_percent/100)*(self.tf-self.t0)

        # if trying to scroll too far right
        elif tf_new >= self.tmax:
            t_shift = direction*(self.tmax - self.tf)

        self.t0 = self.t0 + t_shift
        self.tf = self.tf + t_shift
        self.refresh()

#%% PLUGIN HANDLERS
    
    def plugins_on_traceplotframe_init(self):
        """
        This runs when TracePlotFrame is initialized

        """
        debug_print('\tTracePlotFrame')
        for el in self.plugins:
            try:
                sys.modules['.'.join([el, el])].on_traceplotframe_init(self)
            except:
                if 'on_traceplotframe_init' in dir(sys.modules['.'.join([el, el])]):
                    raise ValueError("Error with on_traceplotframe_init() in " + str(el))
                else:
                    warnings.warn(str(el) + " has no on_traceplotframe_init() method", UserWarning)
        return
    
    def plugins_on_refresh_st_beginning(self):
        """
        This runs at the beginning of refresh_st()

        """
        debug_print('\tTracePlotFrame')
        for el in self.plugins:
            try:
                sys.modules['.'.join([el, el])].on_refresh_st_beginning(self)
            except:
                if 'on_refresh_st_beginning' in dir(sys.modules['.'.join([el, el])]):
                    raise ValueError("Error with on_refresh_st_beginning() in " + str(el))
                else:
                    warnings.warn(str(el) + " has no on_refresh_st_beginning() method", UserWarning)
        return
    
    def plugins_on_refresh_st_end(self):
        """
        This runs at the end of refresh_st()

        """
        debug_print('\tTracePlotFrame')
        for el in self.plugins:
            try:
                sys.modules['.'.join([el, el])].on_refresh_st_end(self)
            except:
                if 'on_refresh_st_end' in dir(sys.modules['.'.join([el, el])]):
                    raise ValueError("Error with on_refresh_st_end() in " + str(el))
                else:
                    warnings.warn(str(el) + " has no on_refresh_st_end() method", UserWarning)
        return
    
    def plugins_on_refresh_beginning(self):
        """
        This runs at the beginning of refresh()

        """
        debug_print('\tTracePlotFrame')
        for el in self.plugins:
            try:
                sys.modules['.'.join([el, el])].on_refresh_beginning(self)
            except:
                if 'on_refresh_beginning' in dir(sys.modules['.'.join([el, el])]):
                    raise ValueError("Error with on_refresh_beginning() in " + str(el))
                else:
                    warnings.warn(str(el) + " has no on_refresh_beginning() method", UserWarning)
        return
    
    def plugins_on_refresh_end(self):
        """
        This runs at the end of refresh()

        """
        debug_print('\tTracePlotFrame')
        for el in self.plugins:
            try:
                sys.modules['.'.join([el, el])].on_refresh_end(self)
            except:
                if 'on_refresh_end' in dir(sys.modules['.'.join([el, el])]):
                    raise ValueError("Error with on_refresh_end() in " + str(el))
                else:
                    warnings.warn(str(el) + " has no on_refresh_end() method", UserWarning)
        return
    
    
#%% MISCELLANEOUS FUNCTIONS

    def fill_fields(self):
        """
        There is a weird bug that deletes info in all the tk entries after about 7 clicks in the trace selection box. 
        I have no idea what causes it. The fields are deleted in the refresh() function on the line that has "self.ax = self.fig.add_subplot(1, 1, 1)".
        It may be an issue with using TkAgg as a backend.
        """
        debug_print('TracePlotFrame')
        # print(self.pref_plot_picks)
        if self.trace_filter_id_entry.get() == '':
            self.detrend1_entry.delete(0, 'end')
            self.detrend2_entry.delete(0, 'end')
            self.f1_entry.delete(0, 'end')
            self.f2_entry.delete(0, 'end')
            self.f3_entry.delete(0, 'end')
            self.trace_filter_id_entry.delete(0, 'end')
            self.trace_filter_dist1_entry.delete(0, 'end')
            self.trace_filter_dist2_entry.delete(0, 'end')
            self.trace_snr_filter_entry.delete(0, 'end')
            
            self.detrend_box.set(self.pref_detrend_type)
            self.detrend1_entry.insert(0, self.pref_detrend_order)
            self.detrend2_entry.insert(0, self.pref_detrend_dspline)
            self.filter_type_box.set(self.pref_filter_type)
            self.f1_entry.insert(0, self.pref_filter_f1)
            self.f2_entry.insert(0, self.pref_filter_f2)
            self.f3_entry.insert(0, self.pref_filter_f3)
            # self.pref_filter_zerophase = False
            self.trace_filter_id_entry.insert(0, self.pref_id_filter_str)
            self.trace_filter_dist1_entry.insert(0, self.pref_dist1_filter)
            self.trace_filter_dist2_entry.insert(0, self.pref_dist2_filter)
            self.trace_snr_filter_entry.insert(0, self.pref_trace_snr)
            
            self.plot_picks_box.set(self.pref_plot_picks)



    def save_figure(self):
        debug_print('TracePlotFrame')
        Path('figures/').mkdir(exist_ok=True)
        fname = 'figures/'+'.'.join(self.fig_name_list)+'.pdf'
        self.fig.savefig(fname)
        
    def get_stream_level_properties(self):
        """
        Gets all GUI fields and stores them in the object

        """
        debug_print('TracePlotFrame')
        # get fields
        # self.pref_plot_picks = self.plot_picks_box.get()
        self.pref_detrend_type = self.detrend_box.get()
        # self.pref_filter_type = self.filter_type_box.get()
        # self.f1_value = float(self.f1_entry.get())
        # self.f2_value = float(self.f2_entry.get())
        # self.f3_value = int(self.f3_entry.get())
        # self.pref_filter_zerophase = bool(self.zerophase_checkbutton_value.get())
        
        # get values from other frames
        # self.efs_parent_dir_new = self.controller.frames[StartPage].efs_parent_dir
        self.efs_parent_dir = self.controller.frames[StartPage].efs_parent_dir
        
        # quick processing
        # if self.pref_detrend_type=='None':
        #     self.pref_detrend = 0
        # elif self.pref_detrend_type=='Demean':
        #     self.pref_detrend = 1
        # elif self.pref_detrend_type=='Detrend':
        #     self.pref_detrend = 2
        return
    
    def get_trace_level_properties(self):
        debug_print('TracePlotFrame')
        return
    
    def print_figure_properties(self):
        """
        This is a function used for debugging figure size issues.

        """
        debug_print('TracePlotFrame')
        bbox = self.ax.get_window_extent().transformed(self.fig.dpi_scale_trans.inverted())
        width, height = bbox.width, bbox.height
        print('----------------------------')
        print('Figure size (in): ', self.fig.get_size_inches())
        print('Figure size (px): ', self.fig.get_size_inches() * self.fig.dpi)
        print('Figure dpi: ', self.fig.dpi)
        print('Axes size (px): ', [width, height])
        print('canvas_widget size: ', [self.canvas_widget.winfo_width(), self.canvas_widget.winfo_height()])
        print('canvas_frame size: ', [self.canvas_frame.winfo_width(), self.canvas_frame.winfo_height()])
        print('canvas size: ', self.canvas.get_width_height())
        print('----------------------------')

        # print(dir(self.canvas))
        
        return
    
    def on_zerophase_checkbutton_toggle(self, var):
        """
        This is a workaround from https://stackoverflow.com/questions/50276202/tkinter-checkbutton-not-working

        """
        debug_print('TracePlotFrame')
        var.set(not var.get())
        # zp_state = self.zerophase_checkbutton['state']
        # print('value: ', self.zerophase_checkbutton_value.get(), ', state: ', zp_state)
        # print(dir(self.zerophase_checkbutton_value))
        
        self.on_filter_option_change(False)
        # self.refresh_st()
        return
    
        
#%% REFRESH FUNCTIONS
        
    def refresh_st(self):
        """
        Loads the file if necessary, filters, updates some preferences, etc.
        Obspy.Stream level operations

        Calls refresh() at the end.
        
        self.st is the original Obspy stream loaded from file(s)
        self.st_subset is the original Obspy stream minus traces filtered out with filter_trace_list()
        self.st_pref is the st_subset with filtered signals
        
        """
        debug_print('TracePlotFrame')
        global st
        global st_orig
        
        
        # get and store properties
        self.get_stream_level_properties()
        self.plugins_on_refresh_st_beginning()
        # print('zerophase: ', self.pref_filter_zerophase)
        # store a few objects
        self.st_pref = self.st_subset.copy()
        self.ntr = len(self.st_subset)
        
        st_orig = self.st_pref.copy()
        
        # if this is a new directory, handle it and populate file listbox. is this necessary/properly located?
        # if self.efs_parent_dir_new != self.efs_parent_dir:
        #     self.efs_parent_dir = self.efs_parent_dir_new
        #     self.populate_file_listbox()
        
        # update filter options
        # self.pref_filter_type = self.filter_type_box.get()

        ##### apply preferences #####
        # if self.pref_detrend==0: # do nothing
        # if self.pref_detrend_type==:
        #     self.st_pref.detrend('simple')
        # elif self.pref_detrend_type==2:
        #     self.st_pref.detrend('linear')
        
        # DETREND OPTION
        if self.pref_detrend_type != 'None':
            self.st_pref.detrend(self.pref_detrend_type.lower(), **self.pref_detrend_options)
        
        # filter and set filter string
        if self.pref_filter_type=='highpass' or self.pref_filter_type=='lowpass':
            self.st_pref.filter(self.pref_filter_type,freq=self.pref_filter_f1, zerophase=self.pref_filter_zerophase)
            self.filter_string = self.pref_filter_type + ", f="+str(self.pref_filter_f1)+" Hz"
            self.name_filter_string = self.pref_filter_type +str(self.pref_filter_f1)
            
        elif self.pref_filter_type=='bandpass':
            self.st_pref.filter(self.pref_filter_type,freqmin=self.pref_filter_f1,freqmax=self.pref_filter_f2,corners=self.pref_filter_f3, zerophase=self.pref_filter_zerophase)
            self.filter_string = self.pref_filter_type + " (fmin,fmax)=("+str(self.pref_filter_f1) + ","+str(self.pref_filter_f2)+") Hz, corners="+str(self.pref_filter_f3)
            self.name_filter_string = self.pref_filter_type +str(self.pref_filter_f1) + "-"+str(self.pref_filter_f2)+"."+str(self.pref_filter_f3)+'corners'
        
        elif self.pref_filter_type=='none':
            self.filter_string = 'unfiltered'
            self.name_filter_string = 'unfiltered'
        
        self.filter_trace_list_snr()
        st = self.st_pref_subset
        
        self.populate_trace_listbox(st)
        
        self.plugins_on_refresh_st_end()
        
        self.new_file = False

        self.refresh()
    
    def refresh(self):
        """
        Refreshes non-stream specific things like plot features.
        
        Obspy.Trace level operations.

        Returns
        -------
        None.

        """
        debug_print('TracePlotFrame')
        
        # if no traces, do nothing
        if self.ntr == 0: 
            print('no traces')
            return
        
        self.fig.clear()
        # debug_print('TracePlotFrame b1')
        # self.ax = self.fig.add_subplot(111)
        self.ax = self.fig.add_subplot(1, 1, 1)
        # debug_print('TracePlotFrame b2')
        
        # clear the previous figure/axes
        # self.ax.clear()
        
        # store some values
        tr = self.st_pref_subset[self.n_current_trace]
        self.tr_pref = tr.copy()
        self.f_nyquist = self.tr_pref.stats['sampling_rate']/2
        # print('nyquist: ', self.f_nyquist)
        # self.nclick = 0
        
        # make tdif compatible with non-efs files
        if hasattr(tr.stats.pick_data,'tdif'):
            tdif = tr.stats.pick_data['tdif']
        else:
            tdif = 0
        
        # update GUI objects
        self.trace_listbox_label['text'] = str(self.n_current_trace+1)+" of "+str(self.ntr)+" traces"
        self.f1_entry['to'] = self.f_nyquist
        self.f2_entry['to'] = self.f_nyquist
        self.f3_entry['to'] = self.f_nyquist
        
        ### update trace-specific values
        # update preferences
        
        # time array relative to eq origin time
        t = self.tr_pref.times() + tdif
        self.tmin = min(t)
        self.tmax = max(t)
        if self.plot_picks_box.get() != '': self.pref_plot_picks = self.plot_picks_box.get()
        # done updating things
        
        self.plugins_on_refresh_beginning()
        
        # this should be added as a plugin
        pickn,pickt,pickq = [],[],[]
        if self.pref_plot_picks=="Real":
            pickn,pickt,pickq = get_picks(self.tr_pref)
            
            # only plot picks if there are picks
            if len(pickt)>0:
                pickt = pickt + tdif
                # setup color array
                pickc = copy.deepcopy(pickn)
                pickc[pickn=='P']=self.p_c
                pickc[pickn=='S']=self.s_c
                
                # find picks of each type with best quality picks
                bpn,bpt,bpq = best_picks(pickn,pickt,pickq)
                p_arrival_time = bpt[bpn=='P']
                s_arrival_time = bpt[bpn=='S']
                
                if len(p_arrival_time)>0: p_arrival_time = p_arrival_time[0]
                if len(s_arrival_time)>0: s_arrival_time = s_arrival_time[0]
                
                # plot picks
                for ii,pt in enumerate(pickt):
                    
                    if pickq[ii]==0.1: 
                        ls = ':'
                        lb = pickn[ii]+' arrival (est.)'
                    else:
                        ls = '-'
                        lb = pickn[ii]+' arrival'
                    if self.pref_plot_picks=='Real':
                        if pickq[ii]>0.1: self.ax.axvline(pt,c=pickc[ii],linestyle=ls,zorder=3,label=lb)
                    # elif self.pref_plot_picks=='All':
                    #     self.ax.axvline(pt,c=pickc[ii],linestyle=ls,zorder=3,label=lb)
        
        # if auto-zoom couldn't find something to zoom to
        if self.tf==-99999:
            self.t0 = t[0]
            self.tf = t[-1]
            self.perc_time_buffer = 1
        # else:
        #     # this is weird, revise
        #     if self.nclick==0:
        #         print(self.t0,self.tf)
        #         self.nclick = 1
            
        # debug_print('TracePlotFrame g')    
        
        d = self.tr_pref.data
        self.ax.plot(t,d,c='k',zorder=1,linewidth=0.5,label="Trace")
        # self.fig.tight_layout(pad=1.4)
        
        # self.fig.subplots_adjust(left=0.05, bottom=0.07, right=0.95, top=0.95, wspace=0, hspace=0)
        
        #%% determine and set good xlim and ylim based on what is plotted
        delta_t = self.tf-self.t0
        
        # time_buffer is the amount of time on edges of trace to exclude from dmax calculation
        time_buffer = (self.perc_time_buffer/100)*delta_t
        # print('t window: ',self.t0+time_buffer,self.tf-time_buffer)
        dwin = d[np.logical_and((t>=self.t0+time_buffer),(t<=self.tf-time_buffer))]
        
        d_mean = np.nanmean(dwin)
        self.dmax = max(abs(dwin-d_mean))
        
        ylimits = (d_mean + np.array([-1,1],dtype=np.float64) * self.dmax * 1.2)
        self.ax.set_ylim(ylimits)
        
        xlimits = [self.t0,self.tf]
        
        self.ax.set_xlim(xlimits)

        # debug_print('TracePlotFrame h')
        
        #%%
        ## plot stuff
        self.ax.legend(loc='lower right',fontsize=self.infobox_text_size,framealpha=0.9,frameon=True)
        
        self.fig_title_list = [str(el) for el in [tr.stats.event_data['evid'],tr.id,"M"+str(np.round(tr.stats.event_data['qmag1'],2)),self.filter_string]]
        self.title_str = ' | '.join(self.fig_title_list)
        
        
        self.fig_name_list = [str(el) for el in [tr.stats.event_data['evid'],tr.id,self.name_filter_string,str(int(self.t0))+'s-'+str(int(self.tf))+'s']]
        
        
        self.ax.set_title(self.title_str,fontsize=self.title_size)
        props = dict(boxstyle='round', facecolor='white', alpha=0.9)
        self.ax.annotate(tr.id,xy=(0,1),xytext=(5,-5),fontsize=self.infobox_text_size,xycoords='axes fraction', textcoords='offset points',bbox=props,horizontalalignment='left',verticalalignment='top',zorder=2000)
        self.ax.tick_params(axis='both', which='major', labelsize=self.axes_ticklabel_size)
        self.ax.set_xlabel("Time relative to origin (s)", fontsize=self.axes_label_size)
        
        self.ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%3.1e'))
        
        # print(dir(self.ax))
        
        ### upper right info
        if self.pref_plot_picks=="Real":
            # print('creating pick labels')
            pick_textlabels = "\n".join([el.upper()+" arrival: "+str(np.round(pickt[ii],2))+" s" for ii,el in enumerate(pickn)])
            if len(pick_textlabels)>0: pick_textlabels = "\n"+pick_textlabels
        else: 
            pick_textlabels = ''
        ep_dist = tr.stats.station_data.deldist
        textstr = '\n'.join(("Ep dist: "+str(round(ep_dist,2))+" km","Fs: "+str(round(tr.stats.sampling_rate,1))+" Hz"))+pick_textlabels
        props = dict(boxstyle='round', facecolor='white', alpha=0.9)
        self.ax.annotate(textstr,xy=(1,1),xytext=(-5,-5),fontsize=self.infobox_text_size,xycoords='axes fraction', textcoords='offset points',bbox=props,horizontalalignment='right',verticalalignment='top',zorder=2000)
        
        self.plugins_on_refresh_end()
        
        self.update_external_figures()
        
        # this fixes figure size. probably temporary fix but it works for now
        self.fig.set_dpi(100)
        self.fig.tight_layout()
        # pos1 = self.ax.get_position() # get the original position 
        # pos2 = [pos1.x0 + 0.1, pos1.y0,  pos1.width, pos1.height]
        # self.ax.set_position(pos2)
        # print(pos1)
        
        # trmax = np.abs(tr.max())
        # trmean = np.mean(np.abs(tr.data))
        # print("SNR: %6.3f" %(trmax/trmean))
        
        self.update()
        self.canvas.draw()
        self.fill_fields()
        # print(sys.exc_info()[2])
        
        # print(plt.get_fignums())
        # self.update()
        # self.print_figure_properties()

def plugins_bind_events(frame):
    global plugins
    debug_print('TracePlotFrame')
    
    for el in plugins:
        try:
            sys.modules['.'.join([el, el])].on_bind_init(frame)
        except:
            warnings.warn("Error binding events for " + str(el) + ". Skipping.", UserWarning)
    return

#%% Program starts here

#%% LOAD SETTINGS
global config_debug_mode
global plugins
global app_geometry
global active_plugins
global frame_tpf

load_config_handler()
plugins = load_plugins(active_plugins)

app_geometry = np.add(DEF_GEOMETRY, [plugin_total_height, 0])

app = SeisPlot()
frame_tpf = app.frames[TracePlotFrame]
frame_start = app.frames[StartPage]

# bind events to functions
frame_tpf.detrend_box.bind('<<ComboboxSelected>>', lambda _: frame_tpf.on_detrend_option_change())
frame_tpf.filter_type_box.bind('<<ComboboxSelected>>', lambda _: frame_tpf.on_filter_option_change(False))
frame_tpf.f1_entry.bind('<FocusOut>', lambda _: frame_tpf.on_filter_option_change(False))
frame_tpf.f2_entry.bind('<FocusOut>', lambda _: frame_tpf.on_filter_option_change(False))
frame_tpf.f3_entry.bind('<FocusOut>', lambda _: frame_tpf.on_filter_option_change(False))
frame_tpf.plot_picks_box.bind('<<ComboboxSelected>>', lambda _: frame_tpf.refresh())

frame_tpf.file_listbox.bind('<<ListboxSelect>>', lambda _: frame_tpf.on_file_change())
frame_tpf.trace_listbox.bind('<<ListboxSelect>>', lambda _: frame_tpf.on_trace_change())

# trace detrending
frame_tpf.detrend1_entry.bind('<FocusOut>', lambda _: frame_tpf.on_detrend_option_change())
frame_tpf.detrend2_entry.bind('<FocusOut>', lambda _: frame_tpf.on_detrend_option_change())

# trace filtering
frame_tpf.trace_filter_id_entry.bind('<FocusOut>', lambda _: frame_tpf.on_distance_or_id_filter_change())
frame_tpf.trace_filter_dist1_entry.bind('<FocusOut>', lambda _: frame_tpf.on_distance_or_id_filter_change())
frame_tpf.trace_filter_dist2_entry.bind('<FocusOut>', lambda _: frame_tpf.on_distance_or_id_filter_change())
frame_tpf.trace_snr_filter_entry.bind('<FocusOut>', lambda _: frame_tpf.on_snr_filter_change())

# frame_tpf.canvas.bind()
frame_tpf.canvas.callbacks.connect('button_press_event', frame_tpf.on_canvas_click)

plugins_bind_events(frame_tpf)

# stop script on window close
app.protocol("WM_DELETE_WINDOW", app.quit_program)

app.geometry(get_geometry_string(app_geometry))
app.minsize(width=app_geometry[1], height=app_geometry[0])

app.mainloop()
