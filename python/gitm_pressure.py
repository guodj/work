#!/home/guod/anaconda3/bin/python
#-------------------------------------------------------------------------------
#calc_pressure() - calculate the pressure in gitm
#-------------------------------------------------------------------------------
import gitm
import numpy as np
from spacepy.datamodel import dmarray
import gitm_create_coordinate as gcc
import gitm_3D_const_alt as g3ca
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

def calc_pressure(g):
    '''
    Calculate the pressure.
    Input: g          = gitm data
    Output: g with `pressure` added
    '''
    Re = 6371*1000 # Earth radius, unit: m
    k = 1.38064852*1e-23
    n = g['O(!U3!NP)']+g['O!D2!N']+g['N!D2!N']+g['N(!U4!NS)']+g['NO']+g['He']
    g['pressure'] = dmarray(
            n*k*g['Temperature'],
            attrs={'units':'kg/(m*s^2)', 'scale':'log', 'name':'pressure'})
    return g


if __name__ == '__main__':
    None
