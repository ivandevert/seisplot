#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 18 16:29:23 2022

@author: ivandevert
"""
from tkinter import ttk
import numpy as np
from matplotlib import pyplot as plt

def on_grid(obj, frame):
    obj.sample_button = ttk.Button(frame, text="Sample button",command= lambda:print_test('This is a test button and function'))
    obj.sample_button.pack(side='left')
    
    obj.sample_plot_button = ttk.Button(frame, text="Plot button",command= lambda:plot_test(obj))
    obj.sample_plot_button.pack(side='right')
    return

def on_traceplotframe_init(obj):
    
    return

def on_refresh_st_beginning(obj):
    
    return

def on_refresh_st_end(obj):
    
    return

def on_refresh_beginning(obj):
    
    return

def on_refresh_end(obj):
    
    return

def on_bind_init(obj):
    # obj.*tkobject*.bind()
    return

def plot_test(obj):
    x = np.linspace(-100, 1000, 10000)
    obj.ax.plot(x, obj.dmax * 0.8 * np.random.rand() * np.sin((x)/2 + np.random.rand()*2*np.pi))
    
    obj.update()
    obj.canvas.draw()

def print_test(string):
    print(string)
    return