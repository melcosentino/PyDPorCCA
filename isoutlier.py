import numpy as np
import scipy.special as special

K = -1/(np.sqrt(2) * special.erfcinv(3/2))

def isoutliers(x):
    """
    Identifies and locates outliers in the data. Translated from Matlab function with the same name by Clea Parcerisas
    """

    mad = K * np.nanmean(np.abs(x - np.nanmean(x)))
    return np.abs(x - np.nanmean(x)) > 3*mad

