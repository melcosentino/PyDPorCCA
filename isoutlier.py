import numpy as np
import scipy.special as special


K = -1/(np.sqrt(2) * special.erfcinv(3/2))


def isoutliers(x):
    mad = K * np.median(np.abs(x - np.median(x)))
    return np.abs(x - np.median(x)) > 3*mad

