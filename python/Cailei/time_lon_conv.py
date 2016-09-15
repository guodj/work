from datetime import time
from datetime import datetime
from datetime import timedelta

def time_to_lon(t):
    return t.hour * 15 + t.minute * 0.25 + t.second / 240 + t.microsecond / 240000000

def lt_to_lon(ut, lt):
    lon = time_to_lon(lt) - time_to_lon(ut)
    if (lon < 0):
        lon = 360 - lon
    return lon


