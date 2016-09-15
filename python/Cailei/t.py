import numpy as np
import matplotlib.pyplot as plt
from gitm_routines_new import *

header = read_gitm_file_header('F:\\data\\gitm\\3DALL_t000318_000000.bin')
i = 0
for var in header['vars']:
    print(i, header['vars'][i].decode())
    i += 1

vars = [0,1,2,3,15,18]
data = get_gitm_vars_range('F:\\data\\gitm\\3DALL_t000318_000000.bin', vars)

i += 1
