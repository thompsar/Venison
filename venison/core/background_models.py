import numpy as np

def bg3D(x,b,c):
    return b*np.exp(-x/c)

def bgFractal(x,b,c,d):
    return b*np.exp(-(np.sign(x)*np.abs(x)**(np.abs(d)/3))/c)