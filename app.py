import sys

import numpy as np
import panel as pn
import hvplot.pandas
import holoviews as hv
import pandas as pd
import holoviews.operation.datashader as hd
import datashader as ds
from holoviews.streams import RangeXY
import tskit
import utils

# Usage: panel serve app.py --args /path/to/trees-file
path = sys.argv[1]
ts = tskit.load(path)
ti = utils.TreeInfo(ts, 1)


def mutations_data(ts, log_time_offset=0.01):
    mutations = ts.tables.mutations
    mutations_time = ts.mutations_time.copy()
    mutations_node = ts.mutations_node.copy()
    unknown = tskit.is_unknown_time(mutations_time)
    mutations_time[unknown] = ts.nodes_time[mutations_node[unknown]]
    mutations_log_time = np.log10(mutations_time + log_time_offset)
    node_flag = ts.nodes_flags[mutations_node]
    position = ts.sites_position[mutations.site]
    array = np.column_stack(
        [position, mutations_node, mutations_time, mutations_log_time, node_flag]
    )
    df = pd.DataFrame(
        array, columns=["position", "mutation_node", "time", "log_time", "node_flag"]
    )
    df = df.astype(
        {
            "position": "float64",
            "mutation_node": "int",
            "time": "float64",
            "log_time": "float64",
            "mutation_node": "int",
            "node_flag": "int",
        }
    )
    return pd.DataFrame(df)


df_mutations = mutations_data(ts)


def filter_points(points, x_range, y_range):
    if x_range and y_range:
        return points[x_range, y_range]
    return points


def hover_points(points, threshold=5000):
    ### Return points to interact with via hover if there are fewer than threshold
    if len(points) > threshold:
        return points.iloc[:0]
    return points


def shaded_points(points, threshold=5000):
    ### Return points to datashade if there are more than threshold
    if len(points) > threshold:
        return points
    return points.iloc[:0]


def make_hist_on_axis(dimension, points, num_bins=30):
    ### Make histogram function for a specified axis of a scatter plot
    def compute_hist(x_range, y_range):
        filtered_points = filter_points(points, x_range, y_range)
        hist = hv.operation.histogram(
            filtered_points, dimension=dimension, num_bins=num_bins, normed="height"
        )
        return hist

    return compute_hist


def make_hist(data, title, bins_range, log_y=True, plot_width=800):
    ### Make histogram from given count data
    count, bins = np.histogram(data, bins=bins_range)
    ylabel = "log(Count)" if log_y else "Count"
    np.seterr(divide="ignore")
    if log_y:
        count = np.log10(count)
        count[count == -np.inf] = 0
    histogram = hv.Histogram((count, bins)).opts(
        title=title, ylabel=ylabel, tools=["hover"]
    )
    histogram = histogram.opts(shared_axes=False, width=round(plot_width / 2))
    return histogram


def make_hist_panel(log_y, plot_width=800):
    ### Make row of mhistograms for holoviews panel
    overall_site_hist = make_hist(
        ti.sites_num_mutations,
        "Mutations per site",
        range(29),
        log_y=log_y,
        plot_width=plot_width,
    )
    overall_node_hist = make_hist(
        ti.nodes_num_mutations,
        "Mutations per node",
        range(10),
        log_y=log_y,
        plot_width=plot_width,
    )
    return pn.Row(overall_site_hist, overall_node_hist)


def page1():
    return pn.pane.HTML(ts)
    # hv_layout


def page2():
    plot_width = 1000
    log_y_checkbox = pn.widgets.Checkbox(
        name="Log y-axis of Mutations per site/node plots", value=False
    )
    threshold = pn.widgets.IntSlider(
        name="Maximum number of points", start=1000, end=10000, step=100
    )

    points = df_mutations.hvplot.scatter(
        x="position",
        y="time",
        hover_cols=["position", "time", "mutation_node", "node_flag"],
    ).opts(width=plot_width, height=round(plot_width / 2))

    range_stream = hv.streams.RangeXY(source=points)
    streams = [range_stream]

    filtered = points.apply(filter_points, streams=streams)
    time_hist = hv.DynamicMap(
        make_hist_on_axis(dimension="time", points=points, num_bins=10), streams=streams
    )
    site_hist = hv.DynamicMap(
        make_hist_on_axis(dimension="position", points=points, num_bins=10),
        streams=streams,
    )
    hover = filtered.apply(hover_points, threshold=threshold)
    shaded = hd.datashade(filtered, width=400, height=400, streams=streams)

    main = (shaded * hover).opts(
        hv.opts.Points(tools=["hover"], alpha=0.1, hover_alpha=0.2, size=10)
    )

    hist_panel = pn.bind(make_hist_panel, log_y=log_y_checkbox, plot_width=plot_width)

    plot_options = pn.Column(
        pn.pane.Markdown("## Plot Options"),
        log_y_checkbox,
        threshold,
    )

    return pn.Column(main << time_hist << site_hist, hist_panel, plot_options)


pn.extension(sizing_mode="stretch_width")

pages = {"Overview": page1, "Mutations": page2}


def show(page):
    return pages[page]()


starting_page = pn.state.session_args.get("page", [b"Overview"])[0].decode()
page = pn.widgets.RadioButtonGroup(
    value=starting_page,
    options=list(pages.keys()),
    name="Page",
    # sizing_mode="fixed",
    button_type="success",
)
ishow = pn.bind(show, page=page)
pn.state.location.sync(page, {"value": "page"})

ACCENT_COLOR = "#0072B5"
DEFAULT_PARAMS = {
    "site": "Panel Multi Page App",
    "accent_base_color": ACCENT_COLOR,
    "header_background": ACCENT_COLOR,
}
pn.template.FastListTemplate(
    title="As Single Page App",
    sidebar=[page],
    main=[ishow],
    **DEFAULT_PARAMS,
).servable()
