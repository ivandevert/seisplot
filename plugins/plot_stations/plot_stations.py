#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
currently can't resize because of blit





"""

from tkinter import ttk
import numpy as np
from matplotlib import pyplot as plt
import os

from pyrsistent import s

try:
    import cartopy
    import cartopy.io.img_tiles as cimgt
    import cartopy.crs as ccrs
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
except:
    raise ImportError("Package 'cartopy' not installed. Please install.")

# height of the plugin's frame. Suggested: at least 45 per row of buttons, text, etc.
frame_height = 45

plot_stations_fignum = 50

script_path = os.path.dirname(os.path.realpath(__file__))
cachedir = script_path + '/tiles_cache/'

def on_grid(obj, frame):
    # obj.sample_button = ttk.Button(frame, text="Sample button",command= lambda:print_test('This is a test button and function'))
    # obj.sample_button.pack(side='left')
    
    obj.sample_plot_button = ttk.Button(frame, text="Plot stations",command= lambda:on_plot_stations_button_press(obj))
    obj.sample_plot_button.pack(side='left')
    return

def on_traceplotframe_init(obj):
    
    return

def on_refresh_st_beginning(obj):
    
    return

def on_refresh_st_end(obj):
    global annot
    try:
        if hasattr(obj, 'station_plot_fig'):
            obj.PLST_lat, obj.PLST_lon, obj.PLST_ev_lat, obj.PLST_ev_lon = get_latlon(obj.st_pref_subset)

            all_lat, all_lon, junk1, junk2 = get_latlon(obj.st)

            obj.station_plot_fig.clf()

            projection=ccrs.PlateCarree()
            GT = cimgt.Stamen(style='terrain', cache=cachedir)
            # extent = np.array([min(obj.PLST_lon), max(obj.PLST_lon), min(obj.PLST_lat), max(obj.PLST_lat)], dtype=float)
            extent = np.array([min(all_lon), max(all_lon), min(all_lat), max(all_lat)], dtype=float)
            
            dlon = abs(extent[1] - extent[0])
            dlat = abs(extent[3] - extent[2])

            buffer = 0.1 * (dlon + dlat) / 2 * np.array([-1, 1, -1, 1], dtype=float)
            extent = extent + buffer

            obj.PLST_ax = obj.station_plot_fig.add_subplot(111, projection=projection)
            obj.PLST_ax.set_extent(extent)
            
            gl = obj.PLST_ax.gridlines(draw_labels=True, alpha=0.2)
            # gl.xlabels_top = gl.ylabels_right = False
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER

            obj.PLST_img = obj.PLST_ax.add_image(GT, 10, interpolation='spline36', alpha=0.5)

            # does this need to be here?
            # obj.PLST_scatter1, = obj.PLST_ax.plot([], markersize=8, markerfacecolor='r', ls='', mew=0, marker='v')
            # obj.PLST_scatter2, = obj.PLST_ax.plot([], markersize=20, markerfacecolor='k', ls='', mew=0, marker='*')
            # obj.PLST_scatter3, = obj.PLST_ax.plot([], markersize=12, markerfacecolor='b', ls='', mew=0, marker='v')

            obj.PLST_scatter1 = obj.PLST_ax.scatter([]*len(obj.st_pref_subset), []*len(obj.st_pref_subset), s=45, c='r', marker='v')
            obj.PLST_scatter2 = obj.PLST_ax.scatter([], [], s=100, c='k', marker='*')
            obj.PLST_scatter3 = obj.PLST_ax.scatter([], [], s=60, c='b', marker='v')

            annot = obj.PLST_ax.annotate("", xy=(0,0), xytext=(4,4),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"))
            annot.set_visible(False)

            obj.station_plot_fig.canvas.draw()
            
            # cache image
            obj.axbackground = obj.station_plot_fig.canvas.copy_from_bbox(obj.PLST_ax.bbox)
            plt.show(block=False)
    except Exception as err:
        print(err)
        
    return

def on_refresh_beginning(obj):
    
    return

def on_refresh_end(obj):
    if hasattr(obj, 'station_plot_fig'):
        try:
            update_map(obj)
        except Exception as err:
            print(err)

    return

def on_bind_init(obj):
    # obj.*tkobject*.bind()
    return

def on_plot_stations_button_press(obj):
    global TPF
    TPF = obj
    obj.station_plot_fig = plt.figure(plot_stations_fignum) #, figsize=(5,5))
    projection=ccrs.PlateCarree()

    # obj.station_plot_fig, obj.station_plot_ax = plt.subplots(figsize=(10,10), subplot_kw=dict(projection=projection))

    # obj.station_plot_fig = plot_cartopy(lon, lat, 10, 'r', marker='v', projection='local', show=True)

    # obj.station_plot_fig.gca().scatter(ev_lon, ev_lat, 10, 'g', marker='*')
    obj.station_plot_fig.canvas.mpl_connect("motion_notify_event", on_hover)
    obj.refresh_st()
    
    return

def update_map(obj):
    # obj.PLST_scatter1.set_data(obj.PLST_lon, obj.PLST_lat)
    # obj.PLST_scatter2.set_data(obj.PLST_ev_lon, obj.PLST_ev_lat)
    # obj.PLST_scatter3.set_data(obj.tr_pref.stats.station_data['slon'], obj.tr_pref.stats.station_data['slat'])
    offsets = np.hstack((np.array(obj.PLST_lon).reshape(len(obj.PLST_lon), 1), np.array(obj.PLST_lat).reshape(len(obj.PLST_lat), 1)))
    obj.PLST_scatter1.set_offsets(offsets)
    obj.PLST_scatter2.set_offsets([obj.PLST_ev_lon, obj.PLST_ev_lat])
    obj.PLST_scatter3.set_offsets([obj.tr_pref.stats.station_data['slon'], obj.tr_pref.stats.station_data['slat']])

    obj.station_plot_fig.canvas.restore_region(obj.axbackground)

    # obj.PLST_ax.draw_artist(obj.PLST_img)
    obj.PLST_ax.draw_artist(obj.PLST_scatter1)
    obj.PLST_ax.draw_artist(obj.PLST_scatter2)
    obj.PLST_ax.draw_artist(obj.PLST_scatter3)

    obj.station_plot_fig.canvas.blit(obj.PLST_ax.bbox)

    obj.station_plot_fig.canvas.flush_events()
    return

def on_hover(event):
    global annot 
    global TPF
    vis = annot.get_visible()
    if event.inaxes == TPF.PLST_ax:
        cont, ind = TPF.PLST_scatter1.contains(event)
        if cont:
            update_annot(ind, annot)
            annot.set_visible(True)
            TPF.station_plot_fig.canvas.draw_idle()
        else:
            if vis:
                annot.set_visible(False)
                TPF.station_plot_fig.canvas.draw_idle()

    return

def update_annot(ind, annot):
    global TPF
    pos = TPF.PLST_scatter1.get_offsets()[ind["ind"][0]]
    annot.xy = pos
    tr_ind = list(map(str, ind["ind"]))
    
    text = ", ".join(list(set([TPF.st_pref_subset[int(el)].stats["station"] for el in tr_ind])))
    annot.set_text(text)
    # annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
    annot.get_bbox_patch().set_alpha(0.4)

def get_latlon(st):
    lat = np.zeros(len(st))
    lon = np.zeros(len(st))

    for ii, tr in enumerate(st):
        lat[ii] = tr.stats.station_data['slat']
        lon[ii] = tr.stats.station_data['slon']
    
    ev_lat = st[0].stats.event_data['qlat']
    ev_lon = st[0].stats.event_data['qlon']

    return lat, lon, ev_lat, ev_lon

def print_test(string):
    print(string)
    return