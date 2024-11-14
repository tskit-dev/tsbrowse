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


def make_hist_on_axis(dimension, points, x_range=None):
    """
    Make histogram function for a specified axis of a scatter plot
    """

    def compute_hist(x_range, y_range):
        filtered_points = filter_points(points, x_range, y_range)
        opts = {} if x_range is None else {"xlim": x_range}
        hist = hv.operation.histogram(
            filtered_points,
            dimension=dimension,
            bins="auto",
            normed="height",
        ).opts(**opts)
        return hist

    return compute_hist


def make_hist(
    data,
    title,
    bins_range,
    xlabel,
    ylabel="count",
    log_y=True,
    plot_width=None,
    plot_height=None,
):
    """
    Make histogram from given count data
    """
    count, bins = np.histogram(data, bins=bins_range)
    np.seterr(divide="ignore")
    if log_y:
        count = np.log10(count)
        count[count == -np.inf] = 0
        ylabel = f"log({ylabel})"

    histogram = hv.Histogram((count, bins)).opts(
        tools=["hover"],
        shared_axes=False,
        width=plot_width,
        height=plot_height,
        toolbar=None,
        default_tools=[],
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
    )
    return histogram


def filter_hist_data(bp_df, x_range):
    if x_range:
        return bp_df[(bp_df.position >= x_range[0]) & (bp_df.position <= x_range[1])]
    return bp_df


def selected_hist(bp_df):
    def compute_hist(x_range, y_range):
        filtered_hist_data = filter_hist_data(bp_df, x_range)
        trees_hist = filtered_hist_data.hvplot.hist(
            y="position",
            bins="auto",
            normed="height",
        )
        return trees_hist

    return compute_hist


def customise_ticks(plot, element):
    p = plot.state
    first, last = int(np.round(p.y_range.start)), int(np.round(p.y_range.end))
    p.yaxis.ticker = [first, last]
    p.yaxis.major_label_overrides = {first: str(first), last: str(last)}


def center_plot_title(plot, element):
    plot.state.title.align = "center"


def parse_range(range_str):
    try:
        start, end = map(float, range_str.split(":"))
        return (start, end)
    except (ValueError, TypeError):
        return None
