import numpy as np

def bg3D(x,b,c):
    return b*np.exp(-x/c)