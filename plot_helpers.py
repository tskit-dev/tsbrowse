import holoviews as hv
import numpy as np

import config


def hover_points(points, threshold=config.THRESHOLD):
    """
    Return points to interact with via hover if there are fewer than threshold
    """
    if len(points) > threshold:
        return points.iloc[:0]
    return points


def filter_points(points, x_range, y_range):
    if x_range and y_range:
        return points[x_range, y_range]
    return points


def make_hist_matplotlib(data, title, num_bins, log_y=True, xlim=(None, None)):
    """
    Make histogram from given count data using parameters
    suitable for the matplotlib backend
    """
    if xlim[1] is not None:
        data = data[data < xlim[1]]
    if xlim[0] is not None:
        data = data[data > xlim[0]]
    count, bins = np.histogram(data, bins=num_bins)
    ylabel = "log(Count)" if log_y else "Count"
    np.seterr(divide="ignore")
    if log_y:
        count = np.log10(count)
        count[count == -np.inf] = 0

    return hv.Histogram((count, bins)).opts(title=title, ylabel=ylabel)
