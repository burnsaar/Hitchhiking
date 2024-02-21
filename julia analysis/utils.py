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
                   legend = True,
                   edgecolor='k',
                   legend_kwds={'label':'Zoned PDR-2'}
                   )
    plt.legend()
    return fig, ax

def plot_depots(fig, ax, gdf, marker='H', color='aqua', legend_kwds={'label':'Depot'}):
    gdf.plot(ax=ax,
             legend=True,
             markersize=500,
             color=color,
             edgecolor='k',
             marker=marker,
             legend_kwds=legend_kwds)
    ax.legend(legend_kwds)
    return fig, ax