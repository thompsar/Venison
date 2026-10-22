import numpy as np

def bg3D(x,b,c):
    return b*np.exp(-x/c)

# Fit bounds on the Fractal dimensionality d (d = 3 is the 3D model).
FRACTAL_D_BOUNDS = (1.0, 5.0)
FRACTAL_FIT_BOUNDS = [[0, 0, FRACTAL_D_BOUNDS[0]], [np.inf, np.inf, FRACTAL_D_BOUNDS[1]]]

def fractal_d_at_limit(d, tolerance=0.1):
    """True when a fitted d sits at (or within `tolerance` of) a fit bound,
    i.e. the Fractal model is pinned rather than fitting."""
    return min(abs(d - bound) for bound in FRACTAL_D_BOUNDS) < tolerance

# For messages, e.g. "its limit (1–5)".
FRACTAL_D_RANGE_TEXT = "{:g}–{:g}".format(*FRACTAL_D_BOUNDS)

def bgFractal(x,b,c,d):
    return b*np.exp(-(np.sign(x)*np.abs(x)**(np.abs(d)/3))/c)