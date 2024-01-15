import holoviews as hv
import numpy as np
import panel as pn

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


def make_hist_on_axis(dimension, points, num_bins=30):
    """
    Make histogram function for a specified axis of a scatter plot
    """

    def compute_hist(x_range, y_range):
        filtered_points = filter_points(points, x_range, y_range)
        hist = hv.operation.histogram(
            filtered_points, dimension=dimension, num_bins=num_bins, normed="height"
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


def make_hist_panel(tsm, log_y):
    """
    Make row of histograms for holoviews panel
    """
    overall_site_hist = make_hist(
        tsm.sites_num_mutations,
        "Mutations per site",
        range(29),
        log_y=log_y,
        plot_width=config.PLOT_WIDTH,
    )
    overall_node_hist = make_hist(
        tsm.nodes_num_mutations,
        "Mutations per node",
        range(10),
        log_y=log_y,
        plot_width=config.PLOT_WIDTH,
    )
    return pn.Row(overall_site_hist, overall_node_hist)
