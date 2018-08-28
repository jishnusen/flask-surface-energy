import numpy as np


def r_sq(x, y, opt_y):
    residuals = y - opt_y
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    return 1 - (ss_res / ss_tot)
