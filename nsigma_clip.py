from numpy import median,abs

def nsigma_clip(data,nsigma=3.):
    """
      return a boolean array same length as data
      True if data are within nsigma robust standard
        deviations of the median
    """
    d0 = median(data)
    if (len(data)>1): dd0 = 1.48*median(abs(data-d0))
    else: dd0 = 1.

    return abs(data-d0) < nsigma*dd0
