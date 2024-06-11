# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 09:29:49 2024

@author: Aaron
"""


import matplotlib.pyplot as plt

def init_fig(figsize):
    fig, ax = plt.subplots(figsize = figsize)
    plt.xlim([-122.52480, -122.33907])
    plt.ylim([37.70, 37.84])
    plt.xlabel('longitude')
    plt.ylabel('latitude')
    plt.legend()
    return fig, ax

def plot_boundary(fig, ax, SF_boundary):
    SF_boundary.boundary.plot(ax=ax, edgecolor = 'k')    
    return fig, ax

def plot_zoning(fig, ax, SF_zoning):
    SF_zoning.plot(ax=ax, 
                   legend=True,
                   edgecolor='k',
                   hatch='/'
                   ) #                   legend_kwds={'label':'Zoned PDR-2'}
    #ax.legend(['Zoned PDR-2'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize=20) #['Zoned PDR-2'], 
    handles, previous_labels = ax.get_legend_handles_labels()
    return fig, ax, handles, previous_labels

def plot_depots(fig, ax, gdf, handles, previous_labels, marker='p', color='darkorange', legend_kwds={'label':'Depot'}):
    gdf.plot(ax=ax,
             legend=True,
             markersize=1500,
             color=color,
             edgecolor='k',
             marker=marker,
             legend_kwds=legend_kwds)
    
    # ax.legend(['Zoned PDR-2', 'Depot Location'], loc='center left', bbox_to_anchor=(1, 0.5), fontsize=20)
    plt.title('Depot Location and Zoning in San Francisco', fontsize=24)
    
    return fig, ax




    # handles, previous_labels = ax.get_legend_handles_labels() #https://stackoverflow.com/questions/23037548/change-main-plot-legend-label-text
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # if mylabels != None:
    #     plt.legend(title = 'Scenarios:', handles=handles, labels=mylabels, loc='center left', bbox_to_anchor=(0.25, -0.4))
    # plt.figure()