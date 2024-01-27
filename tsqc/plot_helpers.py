import holoviews as hv
import numpy as np

from . import config


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


def make_hist_on_axis(dimension, points):
    """
    Make histogram function for a specified axis of a scatter plot
    """

    def compute_hist(x_range, y_range):
        filtered_points = filter_points(points, x_range, y_range)
        hist = hv.operation.histogram(
            filtered_points, dimension=dimension, bins="auto", normed="height"
        )
        return hist

    return compute_hist


def make_hist(data, title, bins_range, log_y=True, plot_width=800):
    """
    Make histogram from given count data
    """
    count, bins = np.histogram(data, bins=bins_range)
    ylabel = "log(Count)" if log_y else "Count"
    np.seterr(divide="ignore")
    if log_y:
        count = np.log10(count)
        count[count == -np.inf] = 0
    histogram = hv.Histogram((count, bins)).opts(
        title=title, ylabel=ylabel, tools=["hover"]
    )
    histogram = histogram.opts(
        shared_axes=False, width=round(plot_width / 2), toolbar=None, default_tools=[]
    )
    return histogram
