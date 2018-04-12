#!/home/guod/anaconda3/bin/python
#-------------------------------------------------------------------------------
#calc_pressure() - calculate the pressure in gitm
#-------------------------------------------------------------------------------

def calc_pressure(g):
    '''
    Calculate the pressure.
    Input: g          = gitm data
    Output: g with `pressure` added
    '''
    Re = 6371*1000 # Earth radius, unit: m
    k = 1.38064852*1e-23
    n = g['O(!U3!NP)']+g['O!D2!N']+g['N!D2!N']+g['N(!U4!NS)']+g['NO']+g['He']+\
        g['N(!U2!ND)']+g['N(!U2!NP)']+g['H']+g['CO!D2!N']+g['O(!U1!ND)']
    g['pressure'] = n*k*g['Temperature']
    return
