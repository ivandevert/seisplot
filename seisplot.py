#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 08:09:03 2022

@author: ivandevert

something is wrong with picking (maybe in insert_pick_times)

channel naming convention: https://ds.iris.edu/ds/nodes/dmc/data/formats/seed-channel-naming/
    letter order: 
        band code (sampling rate and response band of the instrument)
        instrument code (high/low gain, gravimeter, accelerometer, etc.)
        orientation code (Z/N/E, A/B/C, 1/2/3)
    obspy stream.select docs:
        https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.select.html

TO DO:
    add other filters
    add resampling option
    add error logging
    

"""
import tkinter as tk
from tkinter import ttk
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

efspy_module_parent_dir = "/Users/ivandevert/prog/efspy/resources/"

import sys
sys.path.insert(1,efspy_module_parent_dir)
sys.path.append("..")
sys.path.append("../..")
import numpy as np
# import struct


import obspy
# from obspy import UTCDateTime, geodetics
from resources.EFSpy_module import EFS
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
# import warnings
import copy
from tkinter import filedialog


from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg #, NavigationToolbar2Tk
from matplotlib.figure import Figure

from scipy import signal
# from matplotlib.backend_bases import key_press_handler

LARGE_FONT= ("Verdana", 12)

def slash(pth):
    if len(str(pth))>0:
        if str(pth)[-1]!='/':
            return str(pth)+'/'
        else:
            return str(pth)
    else:
        return '/'
    

def get_picks(tr):
    """

    Parameters
    ----------
    tr : obspy.trace or efs.waveforms[n] (dict)
        trace to get picks from.

    Returns
    -------
    pick_names, pick_times (relative to first sample), pick qualities

    """
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
    
    # for efs objects
    elif type(tr)==dict:
        for pt in picklabels:
            if tr[pt+"name"] != '    ':
                picks = np.append(picks,tr[pt])
                pick_names = np.append(pick_names,tr[pt+"name"].strip())
                pick_qualities = np.append(pick_qualities,float(tr[pt+"q"]))
    return pick_names,picks,pick_qualities

def best_picks(pickn,pickt,pickq):
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
    return [str(path)+el for el in os.listdir(path) if el[-len(filetype):]==filetype]

def popup(title,message_str):
    win = tk.Toplevel()
    win.wm_title(title)

    l = tk.Label(win, text=message_str,anchor='w',justify='left',padx=20,pady=20)
    l.grid(row=0, column=0)

    b = ttk.Button(win, text="Okay", command=win.destroy)
    b.grid(row=1, column=0)
    
def popup_response(frame):
    tp = frame.pref_filter_type
    freq = frame.pref_filter_freq
    f1 = frame.pref_filter_bp_fmin
    f2 = frame.pref_filter_bp_fmax
    corners = frame.pref_filter_bp_corners
    f_nyquist = frame.f_nyquist
    
    # print('DEBUG: ',tp,freq,f1,f2,corners,f_nyquist)
    
    if tp == 'highpass' or tp == 'lowpass':
        btype = tp
        F = freq/f_nyquist
    elif tp == 'bandpass':
        btype = 'band'
        F = [f1/f_nyquist, f2/f_nyquist]
        if F[1]>1: F[1]=1.0
    else:
        print('error showing response')
        return
    
    

    b, a = signal.iirfilter(corners,F,btype=btype,ftype='butter')
    w, h = signal.freqz(b,a,fs=f_nyquist*2,worN=512)
    
    # print(F)
    
    # signal.freqz(b,a,plot=lambda w, h:plt.plot(w, np.abs(h)))
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(w, 20* np.log10(np.abs(h)))
    ax.set_title('Butterworth ' + tp + ' frequency response')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('Amplitude [dB]')
    # ax.axis((0.01, f_nyquist*10, -100, 10))
    ax.grid(which='both', axis='both')
    plt.show()

class TraceViewer(tk.Tk):

    def __init__(self, *args, **kwargs):
        global filepath
        global n_current_file
        
        tk.Tk.__init__(self, *args, **kwargs)

        tk.Tk.wm_title(self, "TraceViewer client")
        
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand = True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        self.frames = {}

        for F in (StartPage, TracePlotFrame):

            frame = F(container, self)

            self.frames[F] = frame

            frame.grid(row=0, column=0, sticky="nsew")

        self.show_frame(StartPage)

    def show_frame(self, cont):
        show = True
        if cont==TracePlotFrame:
            if n_current_file >= 0:
                # print('will show frame')
                show = True
            else:
                show = False
                print("Please select a folder with .efs files.")
        if show:
            frame = self.frames[cont]
            if cont==TracePlotFrame: frame.refresh_st()
            frame.tkraise()
    def quit_program(self):
        self.destroy()
        self.quit()

class StartPage(tk.Frame):

    def __init__(self, parent, controller):
        # print('parent n current: ')
        global n_current_file
        self.controller = controller
        
        
        ### set default preferences
        # self.efs_parent_dir = '/Users/ivandevert/seis/sample_data/efs_out/ridgecrest2019_sample/efs/2019/07/'
        self.efs_parent_dir = '/Volumes/ianhdd/projects/ridgecrest2019/efs/2019/07/'
        
        if not os.path.isdir(self.efs_parent_dir): 
            print("not a dir")
            self.efs_parent_dir = "/"
        
        tk.Frame.__init__(self,parent)
        
        main_frame = tk.Frame(self)
        main_frame.pack()
        dir_select_label = tk.Label(main_frame, text="Directory selection", font=LARGE_FONT)
        

        plot_button = ttk.Button(main_frame, text="Plot traces",command=lambda: self.controller.show_frame(TracePlotFrame) )
        
        quit_button = ttk.Button(main_frame, text='Quit', command=lambda: self.controller.quit_program())
        
        label_efs_parent_dir = tk.Label(main_frame, text="EFS parent directory: ", font=LARGE_FONT)
        self.disp_efs_parent_dir = tk.Label(main_frame,text=self.efs_parent_dir)
        button_change_efs_parent_dir = ttk.Button(main_frame,text='Change path',command=lambda: self.change_efs_parent_dir())
        
        
        dir_select_label.grid(row=0,columnspan=2)
        tk.Label(main_frame,text='').grid(row=2)
        label_efs_parent_dir.grid(row=3,columnspan=2)
        self.disp_efs_parent_dir.grid(row=4,columnspan=2)
        plot_button.grid(row=5,column=0,sticky='e')
        button_change_efs_parent_dir.grid(row=5,column=1,sticky='w')
        tk.Label(main_frame,text='').grid(row=6)
        quit_button.grid(row=7,columnspan=2)
                
        self.update_files()
        
    def change_efs_parent_dir(self):
        
        folder_path = filedialog.askdirectory(initialdir=self.efs_parent_dir,title='Select the folder containing EFS files')
        
        if slash(folder_path) != '/':
            self.efs_parent_dir = slash(folder_path)
            self.update_files()
            self.controller.frames[TracePlotFrame].populate_file_listbox()
        
    def update_files(self):
        # print('update_files()')
        global filepath
        global n_current_file
        # print('DEBUG:',cont.frames)
        # update project path
        self.disp_efs_parent_dir['text'] = self.efs_parent_dir
        
        # update project name        
        self.efs_files = get_files(self.efs_parent_dir,'.efs')
        self.efs_files.sort(reverse=False)
        # print(len(self.efs_files),' efs files found')
        
        if len(self.efs_files)>0:
            n_current_file = 0
            filepath = self.efs_files[n_current_file]
        else:
            n_current_file = -1
            filepath = 'null'
        # print('n_current,filepath: ',n_current_file,filepath)
        

class TracePlotFrame(tk.Frame):
    
    def __init__(self,parent,controller):
        global n_current_file
        # print(dir(parent))
        self.controller = controller
        
        self.efs_files = self.controller.frames[StartPage].efs_files
        self.efs_parent_dir = self.controller.frames[StartPage].efs_parent_dir
        # self.n_current_file = self.controller.frames[StartPage].n_current_file
        
        self.NCOL = 14
        # print("TracePlotFrame.__init__()")
        # change freq defaults to be based on sample rate
        
        ### set default preferences
        self.pref_demean = 1                # (0) none, (1) demean, (2) simple detrend
        self.pref_filter_type = 'highpass'  # 'none', 'lowpass', 'highpass', 'bandpass'
        self.pref_filter_freq = 5           # default freq for lowpass and highpass
        self.pref_filter_bp_fmin = 1        # default fmin for bandpass
        self.pref_filter_bp_fmax = 10       # default fmax for bandpass
        self.pref_filter_bp_corners = 4     # default corners for bandpass
        self.pref_id_filter_str = '*,*,*,*'
        self.pref_dist1_filter = 0
        self.pref_dist2_filter = 99999

        # triggers zooming function
        self.manual_t_range = False
        self.nclick = 1
        self.manual_t_range_str = 'Change time range'
        self.scroll_percent = 10
        
        # plotting colors
        self.p_c = 'r'
        self.s_c = 'c'
        
        # convert pref into str
        demean_values = ['None','Demean','Detrend']
        self.pref_demean_str = demean_values[self.pref_demean]
        
        self.pref_plot_picks = 'Real'       # T/F plot picks
        
        self.DEFAULT_ID_FILTER_STR = '*,*,*,*'
        
        box_width = 10
        col_lb1_width = 12
        
        
        tk.Frame.__init__(self,parent)
        
#%% GUI Layout

        main_frame = tk.Frame(self)
        main_frame.grid(row=0,column=0)
        # self.grid_columnconfigure(0,weight=5)
        # self.grid_rowconfigure(0,weight=5)
        
        #%% parent frames
        navbar_frame = tk.Frame(main_frame,padx=5,pady=5,borderwidth=5,relief='groove')
        file_nav_frame = tk.Frame(main_frame,padx=5,pady=5)
        canvas_frame = tk.Frame(main_frame,padx=5,pady=5,borderwidth=2,relief='groove')
        trace_nav_frame = tk.Frame(main_frame,padx=5,pady=5)
        canvas_nav_frame = tk.Frame(main_frame,padx=5,pady=5)
        pref_frame = tk.Frame(main_frame,padx=5,pady=5,borderwidth=1,relief='groove')
        
        # grid parent frames
        navbar_frame.grid(column=0,row=0,columnspan=3,sticky='ew')
        file_nav_frame.grid(column=0,row=1,rowspan=3,sticky='nw')
        canvas_frame.grid(column=1,row=1,sticky='nsew')
        trace_nav_frame.grid(column=2,row=1,rowspan=3,sticky='ne')
        canvas_nav_frame.grid(column=1,row=2,sticky='s')
        pref_frame.grid(column=1,row=3,sticky='s')
        
        # for r in [0,1,2,3]:
        #     main_frame.grid_rowconfigure(r,weight=1)
        # for c in [0,1,2]:
        #     main_frame.grid_columnconfigure(c,weight=1)
        
        # main_frame.grid_columnconfigure(1,weight=5)
        # main_frame.grid_rowconfigure(1,weight=5)
        
        #%% navigation frame
        frame_label = tk.Label(navbar_frame, text="Trace plot",font=LARGE_FONT)
        frame_label.pack(side='left')
        back_button = ttk.Button(navbar_frame, text='Back',command=lambda: controller.show_frame(StartPage))
        back_button.pack(side='left')
        
        save_button = ttk.Button(navbar_frame,text='Save figure',command=lambda: self.save_figure())
        save_button.pack(side='right')
        
        spectrogram_button = ttk.Button(navbar_frame,text='Plot spectrogram',command=lambda: self.plot_spectrogram())
        ppsd_button = ttk.Button(navbar_frame,text='Plot PSD',command=lambda: self.plot_ppsd())
        
        tk.Label(navbar_frame,text=' ',width=10).pack(side='right')
        spectrogram_button.pack(side='right')
        ppsd_button.pack(side='right')
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
        self.file_listbox = tk.Listbox(file_listbox_frame,height=19,exportselection=False)
        self.file_listbox.pack(side='left',fill='both')
        file_scrollbar = tk.Scrollbar(file_listbox_frame)
        file_scrollbar.pack(side='right',fill='both')
        
        self.file_listbox.config(yscrollcommand=file_scrollbar.set)
        file_scrollbar.config(command=self.file_listbox.yview)
        
        self.file_listbox_label = tk.Label(file_nav_frame,text='file counter')
        self.file_listbox_label.grid(row=3,sticky='n')
        ### end file nav frame
        
        #%% canvas frame
        self.fig = Figure(figsize=(10,5), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=canvas_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill='both',expand=True)
        # self.canvas._tkcanvas.pack()
        ### end canvas frame
        
        #%% canvas navigation frame
        self.zoom_toggle_button = ttk.Button(canvas_nav_frame,text=self.manual_t_range_str, command=lambda: self.zoom_toggle())
        self.zoom_toggle_button.pack(side='left')
        self.zoom_status = tk.Label(canvas_nav_frame,text='',fg='red',width=20)
        self.zoom_status.pack(side='left')
        
        # scroll frame
        self.scroll_frame = tk.Frame(canvas_nav_frame)
        self.scroll_frame.pack(side='right')
        self.scroll_label = tk.Label(self.scroll_frame,text='Scroll')
        self.scroll_left_button = ttk.Button(self.scroll_frame,text='<-',command=lambda: self.scroll(-1),width=2)
        self.scroll_right_button = ttk.Button(self.scroll_frame,text='->',command=lambda: self.scroll(1),width=2)
        self.scroll_left_button.pack(side='left')
        self.scroll_label.pack(side='left')
        self.scroll_right_button.pack(side='left')
        # end scroll frame
        ### end canvas navigation frame
        
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
        self.trace_listbox = tk.Listbox(trace_listbox_frame,height=19,exportselection=False)
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
        
        
        component_var = tk.StringVar()
        dist1_var = tk.StringVar()
        dist2_var = tk.StringVar()
        self.trace_filter_id_entry = ttk.Entry(trace_filter_frame,textvariable=component_var,width=8)
        self.trace_filter_dist1_entry = ttk.Entry(trace_filter_frame,textvariable=dist1_var,width=2)
        self.trace_filter_dist2_entry = ttk.Entry(trace_filter_frame,textvariable=dist2_var,width=5)
        
        self.trace_filter_id_help_button = ttk.Button(trace_filter_frame,text='?',width=0.5,command=lambda: self.filter_help_popup())
        
        trace_filter_frame_label.grid(row=0,column=0,columnspan=4,sticky='w')
        trace_filter_id_label.grid(row=1,column=0,sticky='w')
        trace_filter_distance_label.grid(row=2,column=0,sticky='w')
        trace_filter_to_label.grid(row=2,column=2)
        
        self.trace_filter_id_entry.grid(row=1,column=1,columnspan=2,sticky='w')
        self.trace_filter_id_help_button.grid(row=1,column=3,sticky='e')
        self.trace_filter_dist1_entry.grid(row=2,column=1)
        self.trace_filter_dist2_entry.grid(row=2,column=3)
        
        self.trace_filter_id_entry.insert(0,self.pref_id_filter_str)
        self.trace_filter_dist1_entry.insert(0,str(self.pref_dist1_filter))
        self.trace_filter_dist2_entry.insert(0,str(self.pref_dist2_filter))
        # end trace filtering
        
        
        # # sorting MAYBE COMING SOON
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
        
        # demean row row=0
        demean_values = ['None','Demean','Detrend']
        current_pref_demean = tk.StringVar()
        demean_label = ttk.Label(pref_frame,text='Detrend option: ',width=col_lb1_width,anchor='w')
        self.demean_box = ttk.Combobox(pref_frame, textvariable=current_pref_demean,values=demean_values,state='readonly',width=box_width)
        self.demean_box.set(self.pref_demean_str)
        demean_label.grid(column=0,row=0,sticky='w')
        self.demean_box.grid(column=1,row=0,sticky='w',columnspan=6)
        # end demean row
        
        # filter row row=1
        filter_type_values = ['none','highpass','lowpass','bandpass']
        current_pref_filter_type = tk.StringVar()
        filter_type_label = tk.Label(pref_frame,text="Filter type: ",width=col_lb1_width,anchor='w')
        self.filter_type_box = ttk.Combobox(pref_frame,textvariable=current_pref_filter_type,values=filter_type_values,state='readonly',width=box_width)
        self.filter_type_box.set(self.pref_filter_type)
        
        self.f1_label = tk.Label(pref_frame,text='null',width=col_lb1_width,anchor='e')
        self.f2_label = tk.Label(pref_frame,text='null',width=col_lb1_width,anchor='e')
        self.f3_label = tk.Label(pref_frame,text='null',width=col_lb1_width,anchor='e')
        
        self.f1_entry = ttk.Spinbox(pref_frame,width=5,from_=0.01,to=100.0)
        self.f2_entry = ttk.Spinbox(pref_frame,width=5,from_=0.0,to=100.0)
        self.f3_entry = ttk.Spinbox(pref_frame,width=5,from_=0.0,to=100.0)
        self.show_response_button = ttk.Button(pref_frame,text='Response',command=lambda:popup_response(self))
        
        self.f1_entry.insert(0,str(self.pref_filter_freq))
        self.f2_entry.insert(0,str(self.pref_filter_bp_fmax))
        self.f3_entry.insert(0,str(self.pref_filter_bp_corners))
        
        # grid filter row
        filter_type_label.grid(column=0,row=1,sticky='w')
        self.filter_type_box.grid(column=1,row=1,sticky='w')
        self.f1_label.grid(column=2,row=1,sticky='w')
        self.f1_entry.grid(column=3,row=1,sticky='w')
        self.f2_label.grid(column=4,row=1,sticky='w')
        self.f2_entry.grid(column=5,row=1,sticky='w')
        self.f3_label.grid(column=6,row=1,sticky='w')
        self.f3_entry.grid(column=7,row=1,sticky='w')
        self.show_response_button.grid(column=8,row=1,sticky='e')
        # end filter row
        
        ### plot picks row row=2
        plot_picks_values = ['None','Real']
        plot_picks_label = tk.Label(pref_frame,text='Plot picks: ',width=col_lb1_width,anchor='w')
        plot_picks_label.grid(column=0,row=2,sticky='w')
        current_pref_plot_picks = tk.StringVar()
        self.plot_picks_box = ttk.Combobox(pref_frame,textvariable=current_pref_plot_picks,values=plot_picks_values,state='readonly',width=box_width)
        self.plot_picks_box.grid(column=1,row=2,sticky='w')
        self.plot_picks_box.set(self.pref_plot_picks)
        
        self.populate_file_listbox()
        # end plot picks row
        #%% end pref frame
        
        # for ii in np.arange(2,self.NCOL-2):
        #     self.grid_columnconfigure(ii,weight=1)
#%%
        
        self.n_current_trace = 0
        self.t0 = 0
        self.tf = -99999
        
        if n_current_file >= 0:
            self.on_file_change()
            # self.refresh()
    
    def plot_spectrogram(self):
        tr = self.tr_pref.copy()
        # d = tr.data
        # t = tr.times() + tr.stats.pick_data['tdif']
        
        # fig = plt.figure()
        # ax = fig.add_subplot(1,1,1)
        
        tr.spectrogram()
        
        # obspy.imaging.spectrogram.spectrogram(tr.data,tr.stats.sampling_rate,axes=ax)
        # plt.xticks(ticks=np.arange())
        
        return
    
    def plot_ppsd(self):
        tr = self.tr_pref
        d = tr.data
        # t = tr.times()
        
        
        Pxx, freqs = plt.psd(d,Fs=tr.stats.sampling_rate)
        # fig = plt.figure()
        # ax = fig.add_subplot(1,1,1)
        
        # ax.plot(freqs,Pxx)
        
        # plt.show()
        
        return
    
    def filter_traces(self):
        id_str = self.trace_filter_id_entry.get()
        dist1 = float(self.trace_filter_dist1_entry.get())
        dist2 = float(self.trace_filter_dist2_entry.get())
        
        # if nothing changed, return
        if id_str == self.pref_id_filter_str and dist1 == self.pref_dist1_filter and dist2 == self.pref_dist2_filter: 
            return
        else:
            try:
                network_str,station_str,location_str,channel_str = id_str.lower().split(',')
                self.st_subset = self.st.select(network=network_str,station=station_str,location=location_str,channel=channel_str)
                
                for tr in self.st_subset:
                    deldist = tr.stats.station_data['deldist']
                    if deldist < dist1 or deldist > dist2:
                        self.st_subset.remove(tr)
            except Exception as err:
                print("Invalid string entry. Please try again (click ? for help)")
                print(err)
                return
        self.populate_trace_listbox()
        self.n_current_trace = 0
        self.refresh_st()
        
        return
        
    def filter_help_popup(self):
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
    
    # def on_change_trace_order_press(self):
        
    #     if self.trace_sort_ascend:
    #         self.trace_sort_ascend = False
    #         self.trace_order_text = 'descend'
    #     else:
    #         self.trace_sort_ascend = True
    #         self.trace_order_text = 'ascend'
    #     self.trace_sort_order_button['text'] = 'Sorting by ('+self.trace_order_text+'): '
    
    def populate_trace_listbox(self):
        trace_ids = [str(el).split('|')[0].strip() for el in self.st_subset.traces]
        self.ntr = len(trace_ids)
        self.trace_listbox.delete(0,'end')
        
        if self.ntr > 0:
            for el in trace_ids:
                self.trace_listbox.insert('end',el)
        else:
            self.trace_listbox.insert('end','No traces available')
        self.trace_listbox.select_set(0)
        
    
    def populate_file_listbox(self):
        global n_current_file
        # print('populate_file_listbox()')

        n_current_file = 0
        self.efs_files = self.controller.frames[StartPage].efs_files
        
        filenames = [el.split('/')[-1].split('.')[0] for el in self.efs_files]
        self.nfiles = len(filenames)
        self.file_listbox.delete(0,'end')
        # print(filenames)
        for el in filenames:
            self.file_listbox.insert('end',el)
        self.file_listbox.select_set(0)
        self.file_listbox_label['text'] = str(n_current_file+1)+" of "+str(len(self.efs_files))+" files"
        # self.on_file_change()
    
    def scroll(self,direction):
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
        
        
    def increment_trace(self):
        if self.n_current_trace<self.ntr-1: 
            self.n_current_trace += 1
            self.trace_listbox.selection_clear(0,'end')
            self.trace_listbox.select_set(self.n_current_trace)
        self.manual_t_range = False
        self.reset_zoom()
        self.refresh()
        
    def decrement_trace(self):
        if self.n_current_trace>0: 
            self.n_current_trace -= 1
            self.trace_listbox.selection_clear(0,'end')
            self.trace_listbox.select_set(self.n_current_trace)
        self.manual_t_range = False
        self.reset_zoom()
        self.refresh()
        
    def increment_file(self):
        global n_current_file
        # print('DEBUG: ',n_current_file,len(self.efs_files))
        if n_current_file < len(self.efs_files)-1:
            n_current_file += 1
            self.n_current_trace = 0
            self.t0 = 0
            self.tf = -99999
            self.file_listbox.selection_clear(0,'end')
            self.file_listbox.select_set(n_current_file)
        self.on_file_change()
    
    def decrement_file(self):
        global n_current_file
        if n_current_file > 0:
            n_current_file -= 1
            self.n_current_trace = 0
            self.t0 = 0
            self.tf = -99999
            self.file_listbox.selection_clear(0,'end')
            self.file_listbox.select_set(n_current_file)
        self.on_file_change()
    
    def on_trace_selection_change(self):
        nsel = self.trace_listbox.curselection()[0]
        self.n_current_trace = nsel
        self.refresh()
    
    def on_file_change(self):
        # print('on_file_change()')
        global efs_path
        global n_current_file
        
        self.efs_files = self.controller.frames[StartPage].efs_files
        
        n_current_file = self.file_listbox.curselection()[0]
        self.n_current_trace = 0
        
        self.load_file()
        self.filter_traces()
        
        self.file_listbox_label['text'] = str(n_current_file+1)+" of "+str(len(self.efs_files))+" files"
        
        self.refresh_st()
        self.populate_trace_listbox()
    
    def load_file(self):
        global n_current_file
        
        # this will soon be an issue
        self.filepath = self.efs_files[n_current_file]
        # print('ncurr: ',n_current_file)
        self.efs = EFS(self.filepath,np.float32,np.int32)
        self.efs_pref = copy.deepcopy(self.efs)
        
        self.st = self.efs_pref.to_obspy(keep_pkdata=True)
        self.st_subset = self.st.copy()
    
    def save_figure(self):
        Path('figures/').mkdir(exist_ok=True)
        fname = 'figures/'+'.'.join(self.fig_name_list)+'.pdf'
        self.fig.savefig(fname)
        
        
    def refresh_st(self):
        # print('refresh_st()')
        # self.efs_pref = copy.deepcopy(self.efs)
        # self.st_pref = self.st_subset.copy()
        # this placement isn't ideal
        # update plot_picks
        self.pref_plot_picks = self.plot_picks_box.get()
        
        # self.st = self.efs_pref.to_obspy(keep_pkdata=True)
        self.st_pref = self.st_subset.copy()
        self.ntr = len(self.st_subset)
        
        self.efs_parent_dir_new = self.controller.frames[StartPage].efs_parent_dir
        if self.efs_parent_dir_new != self.efs_parent_dir:
            self.efs_parent_dir = self.efs_parent_dir_new
            self.populate_file_listbox()
            
        # maybe only do this if preferences have changed
        
        ##### get current values #####
        # update pref_demean
        pref_demean_str = self.demean_box.get()
        if pref_demean_str=='None':
            self.pref_demean = 0
        elif pref_demean_str=='Demean':
            self.pref_demean = 1
        elif pref_demean_str=='Detrend':
            self.pref_demean = 2
        
        # update filter options
        self.pref_filter_type = self.filter_type_box.get()
        
        # get values, update labels and entry fields
        if self.pref_filter_type=='none':
            self.f1_label['text'] = '----'
            self.f2_label['text'] = '----'
            self.f3_label['text'] = '----'
            self.f1_entry['state'] = 'disabled'
            self.f2_entry['state'] = 'disabled'
            self.f3_entry['state'] = 'disabled'
            self.show_response_button['state'] = 'disabled'
        elif self.pref_filter_type=='highpass' or self.pref_filter_type=='lowpass':
            self.pref_filter_freq = float(self.f1_entry.get())
            
            self.f2_entry['state'] = 'disabled'
            self.f3_entry['state'] = 'disabled'
            self.f1_label['text'] = 'Frequency: '
            self.f2_label['text'] = '----'
            self.f3_label['text'] = '----'
            self.f1_entry['state'] = 'enabled'
            self.f2_entry['state'] = 'disabled'
            self.f3_entry['state'] = 'disabled'
            self.show_response_button['state'] = 'enabled'

        elif self.pref_filter_type=='bandpass':
            self.pref_filter_bp_fmin = float(self.f1_entry.get())
            self.pref_filter_bp_fmax = float(self.f2_entry.get())
            self.pref_filter_bp_corners = int(self.f3_entry.get())
            
            
            self.f1_label['text'] = 'Freq. min: '
            self.f2_label['text'] = 'Freq. max: '
            self.f3_label['text'] = 'Corners: '
            self.f1_entry['state'] = 'enabled'
            self.f2_entry['state'] = 'enabled'
            self.f3_entry['state'] = 'enabled'
            self.show_response_button['state'] = 'enabled'

        ##### apply preferences #####
        # self.st_pref = self.st.copy()
        # if self.pref_demean==0: # do nothing
        if self.pref_demean==1:
            self.st_pref.detrend('simple')
        elif self.pref_demean==2:
            self.st_pref.detrend('linear')
        
        # filter and set filter string
        if self.pref_filter_type=='highpass' or self.pref_filter_type=='lowpass':
            self.st_pref.filter(self.pref_filter_type,freq=self.pref_filter_freq)
            self.filter_string = self.pref_filter_type + ", f="+str(self.pref_filter_freq)+" Hz"
            self.name_filter_string = self.pref_filter_type +str(self.pref_filter_freq)
        elif self.pref_filter_type=='bandpass':
            self.st_pref.filter(self.pref_filter_type,freqmin=self.pref_filter_bp_fmin,freqmax=self.pref_filter_bp_fmax,corners=self.pref_filter_bp_corners)
            self.filter_string = self.pref_filter_type + " (fmin,fmax)=("+str(self.pref_filter_bp_fmin) + ","+str(self.pref_filter_bp_fmax)+") Hz, corners="+str(self.pref_filter_bp_corners)
            self.name_filter_string = self.pref_filter_type +str(self.pref_filter_bp_fmin) + "-"+str(self.pref_filter_bp_fmax)+"."+str(self.pref_filter_bp_corners)+'corners'
        elif self.pref_filter_type=='none':
            self.filter_string = 'unfiltered'
            self.name_filter_string = 'unfiltered'
            
        self.refresh()
    
    def refresh(self):
        
        # print('refresh()')
        self.ax.clear()
        self.trace_listbox_label['text'] = str(self.n_current_trace+1)+" of "+str(self.ntr)+" traces"
                
        if self.ntr == 0: return
        
        tr = self.st_pref[self.n_current_trace]
        tdif = tr.stats.pick_data['tdif']
        tr_pref = tr.copy()
        
        self.tr_pref = tr_pref
        
        self.f_nyquist = tr_pref.stats['sampling_rate']/2
        
        self.f1_entry['to'] = self.f_nyquist
        self.f2_entry['to'] = self.f_nyquist
        self.f3_entry['to'] = self.f_nyquist
        
        ### update trace-specific values
        # update preferences
        self.pref_plot_picks = self.plot_picks_box.get()
        self.nclick = 0
        
        # time array
        t = tr_pref.times() + tdif
        self.tmin = min(t)
        self.tmax = max(t)
        pickn,pickt,pickq = [],[],[]
        if self.pref_plot_picks!='None':
            pickn,pickt,pickq = get_picks(tr_pref)
            
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
                    elif self.pref_plot_picks=='All':
                        self.ax.axvline(pt,c=pickc[ii],linestyle=ls,zorder=3,label=lb)
        
        # if auto-zoom couldn't find something to zoom to
        if self.tf==-99999:
            self.t0 = t[0]
            self.tf = t[-1]
            self.perc_time_buffer = 1
        else:
            # if two limits have just been chosen, set t0,tf here
            if self.nclick==0:
                # print(self.t0,self.tf)
                self.nclick = 1
            
            
        
        d = tr_pref.data
        self.ax.plot(t,d,c='k',zorder=1,linewidth=0.5,label="Trace")
        self.fig.tight_layout(pad=1.4)

        #%% determine and set good xlim and ylim based on what is plotted
        delta_t = self.tf-self.t0
        time_buffer = (self.perc_time_buffer/100)*delta_t
        # print('t window: ',self.t0+time_buffer,self.tf-time_buffer)
        dwin = d[np.logical_and((t>=self.t0+time_buffer),(t<=self.tf-time_buffer))]
        d_mean = np.mean(dwin)
        dmax = max(abs(dwin-d_mean))
        ylimits = (d_mean + np.array([-1,1],dtype=np.float64) * dmax * 1.2)
        self.ax.set_ylim(ylimits)
        
        xlimits = [self.t0,self.tf]
        self.ax.set_xlim(xlimits)
        
        #%%
        ## plot stuff
        self.ax.legend(loc='lower right',fontsize=8,framealpha=0.9,frameon=True)
        
        self.fig_title_list = [str(el) for el in [tr.stats.event_data['evid'],tr.id,"M"+str(np.round(tr.stats.event_data['qmag1'],2)),self.filter_string]]
        self.title_str = ' | '.join(self.fig_title_list)
        
        
        self.fig_name_list = [str(el) for el in [tr.stats.event_data['evid'],tr.id,self.name_filter_string,str(int(self.t0))+'s-'+str(int(self.tf))+'s']]
        
        
        self.ax.set_title(self.title_str,fontsize=10)
        props = dict(boxstyle='round', facecolor='white', alpha=0.9)
        self.ax.annotate(tr.id,xy=(0,1),xytext=(5,-5),fontsize=8,xycoords='axes fraction', textcoords='offset points',bbox=props,horizontalalignment='left',verticalalignment='top',zorder=2000)
        self.ax.tick_params(axis='both', which='major', labelsize=8)
        self.ax.set(xlabel="Time relative to origin (s)")
        
        ### upper right info
        if self.pref_plot_picks!='None':
            pick_textlabels = "\n".join([el.upper()+" arrival: "+str(np.round(pickt[ii],2))+" s" for ii,el in enumerate(pickn)])
            if len(pick_textlabels)>0: pick_textlabels = "\n"+pick_textlabels
        else: 
            pick_textlabels = ''
        ep_dist = tr.stats.station_data.deldist
        textstr = '\n'.join(("Ep dist: "+str(round(ep_dist,2))+" km","Fs: "+str(round(tr.stats.sampling_rate,1))+" Hz"))+pick_textlabels
        props = dict(boxstyle='round', facecolor='white', alpha=0.9)
        self.ax.annotate(textstr,xy=(1,1),xytext=(-5,-5),fontsize=8,xycoords='axes fraction', textcoords='offset points',bbox=props,horizontalalignment='right',verticalalignment='top',zorder=2000)
        
        self.canvas.draw()
        self.update()
        
    def on_click(self, event):
        # print(self.nclick,[event.xdata,event.ydata])
        # print(event)
        # if on manual zoom
        if self.manual_t_range:
            if event.inaxes is not None:
                # print('click')
                if self.nclick==1:
                    self.t1 = event.xdata
                    self.nclick += 1
                    self.zoom_status['text'] = 'Select t2'
                elif self.nclick==2:
                    self.t2 = event.xdata
                    trange = [self.t1,self.t2]
                    trange.sort()
                    # print('trange: ',trange)
                    self.t0 = trange[0]
                    self.tf = trange[1]
                    self.nclick = 0
                    self.zoom_status['text'] = 'Select t1'
                    self.refresh()
                
            else:
                print('clicked outside of axes bounds but inside plot window')
    
    def reset_zoom(self):
        self.tf = -99999
        
    def zoom_toggle(self):
        
        # changing to auto t range
        if self.manual_t_range:
            self.manual_t_range = False
            self.zoom_toggle_button['text'] = 'Change time range'
            self.zoom_status['text'] = ''
            self.reset_zoom()
            
        # changing to manual t range
        else:
            self.manual_t_range = True
            self.zoom_toggle_button['text'] = 'Full time range'
            self.zoom_status['text'] = 'Select t1'
            
        self.refresh()

app = TraceViewer()
frame_tpf = app.frames[TracePlotFrame]
frame_start = app.frames[StartPage]

# bind events to functions
frame_tpf.demean_box.bind('<<ComboboxSelected>>', lambda _: frame_tpf.refresh_st())
frame_tpf.filter_type_box.bind('<<ComboboxSelected>>', lambda _: frame_tpf.refresh_st())
frame_tpf.f1_entry.bind('<FocusOut>', lambda _: frame_tpf.refresh_st())
# frame_tpf.filter_bp_fmin_entry.bind('<FocusOut>', lambda _: frame_tpf.refresh_st())
frame_tpf.f2_entry.bind('<FocusOut>', lambda _: frame_tpf.refresh_st())
frame_tpf.f3_entry.bind('<FocusOut>', lambda _: frame_tpf.refresh_st())
frame_tpf.plot_picks_box.bind('<<ComboboxSelected>>', lambda _: frame_tpf.refresh_st())

frame_tpf.file_listbox.bind('<<ListboxSelect>>', lambda _: frame_tpf.on_file_change())
frame_tpf.trace_listbox.bind('<<ListboxSelect>>', lambda _: frame_tpf.on_trace_selection_change())

# trace filtering
frame_tpf.trace_filter_id_entry.bind('<FocusOut>', lambda _: frame_tpf.filter_traces())
frame_tpf.trace_filter_dist1_entry.bind('<FocusOut>', lambda _: frame_tpf.filter_traces())
frame_tpf.trace_filter_dist2_entry.bind('<FocusOut>', lambda _: frame_tpf.filter_traces())

# frame_tpf.canvas.bind()
frame_tpf.canvas.callbacks.connect('button_press_event', frame_tpf.on_click)


# stop script on window close
app.protocol("WM_DELETE_WINDOW", app.quit_program)

app.mainloop()





