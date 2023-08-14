import sys
import logging

import numpy as np
import panel as pn
import hvplot.pandas
import holoviews as hv
import pandas as pd
import holoviews.operation.datashader as hd
import tskit
import utils
import bokeh.models as bkm
import pathlib
import functools

logger = logging.Logger(__file__)

# Usage: panel serve app.py --args /path/to/trees-file
path = pathlib.Path(sys.argv[1])
trees_file = path.name
logger.warning(f"Loading {path}")
ts = tskit.load(path)
ti = utils.TreeInfo(ts, 1)

# NOTE using "warning" here so that we can get some output
# from them. Will need to do this better at some point,
# with configurable output levels.
logger.warning(f"Computing mutations data frame")
df_mutations = ti.mutations_data()
logger.warning(f"Computing edges data frame")
df_edges = ti.edges_data()
logger.warning(f"Computing Trees data frame")
df_trees = ti.trees_data()

# Global plot settings
plot_width = 1000
plot_height = 600
threshold = 1000  # max number of points to overlay on a plot


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


def make_hist_matplotlib(data, title, num_bins, log_y=True, xlim=(None, None)):
    ### Make histogram from given count data using parameters suitable for the matplotlib backend
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


def make_hist_panel(log_y):
    ### Make row of histograms for holoviews panel
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
    hv.extension("bokeh")
    plot_width = 1000
    log_y_checkbox = pn.widgets.Checkbox(
        name="Log y-axis of Mutations per site/node plots", value=False
    )

    points = df_mutations.hvplot.scatter(
        x="position",
        y="time",
        hover_cols=["position", "time", "mutation_node", "node_flag"],
    ).opts(width=plot_width, height=plot_height)

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

    hist_panel = pn.bind(make_hist_panel, log_y=log_y_checkbox)

    plot_options = pn.Column(
        pn.pane.Markdown("## Plot Options"),
        log_y_checkbox,
    )

    return pn.Column(main << time_hist << site_hist, hist_panel, plot_options)


def page3():
    hv.extension("bokeh")
    df_edges["parent_time_right"] = df_edges["parent_time"]
    lines = hv.Segments(
        df_edges, kdims=["left", "parent_time", "right", "parent_time_right"]
    )
    range_stream = hv.streams.RangeXY(source=lines)
    streams = [range_stream]
    filtered = lines.apply(filter_points, streams=streams)
    hover = filtered.apply(hover_points, threshold=threshold)
    shaded = hd.datashade(filtered, streams=streams)
    hover_tool = bkm.HoverTool(
        tooltips=[
            ("child", "@child"),
            ("parent", "@parent"),
            ("span", "@span"),
            ("branch_length", "@branch_length"),
        ]
    )
    main = (shaded * hover).opts(
        hv.opts.Segments(
            tools=[hover_tool],
            width=plot_width,
            height=plot_height,
            xlabel="Position",
            ylabel="Time",
        )
    )

    return pn.Column(main)


def page4():
    hv.extension("bokeh")
    node_id_input = pn.widgets.TextInput(value="", name="Node ID")
    df_edges["parent_time_right"] = df_edges["parent_time"]
    tabulator = pn.widgets.Tabulator(show_index=False)

    def plot_data(node_id):
        if len(node_id) > 0:
            filtered_df = df_edges[df_edges["child"] == int(node_id)]
            segments = hv.Segments(
                filtered_df,
                kdims=["left", "parent_time", "right", "parent_time_right"],
                vdims=["child", "parent", "span", "branch_length"],
            )
            hover_tool = bkm.HoverTool(
                tooltips=[
                    ("child", "@child"),
                    ("parent", "@parent"),
                    ("span", "@span"),
                    ("branch_length", "@branch_length"),
                ]
            )
            segments = segments.opts(
                width=plot_width,
                height=plot_height,
                tools=[hover_tool],
                xlabel="Position",
                ylabel="Time",
            )

            filtered_df = filtered_df.drop(columns=["parent_time_right"])
            tabulator.value = filtered_df

            return segments
        else:
            return pn.pane.Markdown("Please enter a Node ID.")

    dynamic_plot = pn.bind(plot_data, node_id=node_id_input)

    return pn.Column(node_id_input, dynamic_plot, tabulator)


def page5():
    hv.extension("matplotlib")
    max_num_sites = int(df_trees.num_sites.max())
    min_num_sites = int(df_trees.num_sites.min())
    sites_default_num_bins = min(int(np.sqrt(max_num_sites)), 50)

    spans = df_trees.right - df_trees.left
    min_span = int(np.min(spans))
    max_span = int(np.max(spans))
    spans_default_num_bins = min(50, int(np.sqrt(max_span)))

    sites_hist_func = functools.partial(
        make_hist_matplotlib,
        df_trees.num_sites,
        "Sites per tree",
        num_bins=sites_default_num_bins,
        log_y=True,
        xlim=(min_num_sites, max_num_sites),
    )

    spans_hist_func = functools.partial(
        make_hist_matplotlib,
        spans,
        "Genomic span per tree",
        num_bins=spans_default_num_bins,
        log_y=True,
        xlim=(min_span, max_span),
    )

    if len(df_trees) > 1:
        sites_log_checkbox = pn.widgets.Checkbox(
            name="log y-axis", value=True, align=("center")
        )

        sites_xlim_slider = pn.widgets.IntRangeSlider(
            name="x limits",
            start=min_num_sites,
            value=(min_num_sites, max_num_sites),
            end=max_num_sites,
            step=10,  # TODO set a sensible step size
            align=("center"),
        )

        sites_bins_slider = pn.widgets.IntSlider(
            name="Number of bins",
            start=1,
            end=max(500, int(max_num_sites)),  # TODO set a sensible end value
            value=sites_default_num_bins,
            step=1,  # TODO set a sensible step size
            align=("center"),
        )

        sites_reset_checkbox = pn.widgets.Checkbox(
            name="reset", value=True, align=("center")
        )

        def reset_sites(reset=False):
            if reset:
                sites_log_checkbox.value = True
                sites_xlim_slider.value = (min_num_sites, max_num_sites)
                sites_bins_slider.value = sites_default_num_bins
                sites_reset_checkbox.value = False

        sites_reset_panel = pn.bind(reset_sites, reset=sites_reset_checkbox)

        sites_hist_panel = pn.bind(
            sites_hist_func,
            log_y=sites_log_checkbox,
            xlim=sites_xlim_slider,
            num_bins=sites_bins_slider,
        )

        spans_log_checkbox = pn.widgets.Checkbox(
            name="log y-axis", value=True, align=("center")
        )

        spans_xlim_slider = pn.widgets.IntRangeSlider(
            name="Maximum tree span",
            start=min_span,
            value=(min_span, max_span),
            end=max_span,
            step=1_000,  # TODO set a sensible step size
            align=("center"),
        )

        spans_bins_slider = pn.widgets.IntSlider(
            name="Number of bins",
            start=1,
            end=min(int(max_span), 500),  # TODO set a sensible end value
            value=spans_default_num_bins,
            step=1,  # TODO set a sensible step size
            align=("center"),
        )

        spans_reset_checkbox = pn.widgets.Checkbox(
            name="reset", value=True, align=("center")
        )

        def reset_spans(reset=False):
            if reset:
                spans_log_checkbox.value = True
                spans_xlim_slider.value = (min_span, max_span)
                spans_bins_slider.value = spans_default_num_bins
                spans_reset_checkbox.value = False

        spans_reset_panel = pn.bind(reset_spans, reset=spans_reset_checkbox)

        spans_hist_panel = pn.bind(
            spans_hist_func,
            log_y=spans_log_checkbox,
            xlim=spans_xlim_slider,
            num_bins=spans_bins_slider,
        )

        hist_panel = pn.Row(
            pn.Column(
                sites_reset_panel,
                sites_hist_panel,
                pn.Column(
                    pn.pane.Markdown("Plot Options:", align=("center")),
                    sites_reset_checkbox,
                    sites_log_checkbox,
                    sites_xlim_slider,
                    sites_bins_slider,
                ),
            ),
            pn.Column(
                spans_reset_panel,
                spans_hist_panel,
                pn.Column(
                    pn.pane.Markdown("Plot Options:", align=("center")),
                    spans_reset_checkbox,
                    spans_log_checkbox,
                    spans_xlim_slider,
                    spans_bins_slider,
                ),
            ),
        )
    else:
        sites_hist_single_tree = make_hist_matplotlib(
            df_trees.num_sites, "Sites per tree", num_bins="auto", log_y=False
        )
        sites_hist_single_tree.opts(xticks=[df_trees.num_sites[0]], yticks=[1])
        spans_hist_single_tree = make_hist_matplotlib(
            spans, "Genomic span per tree", num_bins="auto", log_y=False
        )
        spans_hist_single_tree.opts(xticks=[spans[0]], yticks=[1])
        hist_panel = pn.Row(
            pn.Column(sites_hist_single_tree), pn.Column(spans_hist_single_tree)
        )

    return hist_panel


pn.extension(sizing_mode="stretch_width")
pn.extension("tabulator")

pages = {
    "Overview": page1,
    "Mutations": page2,
    "Edges": page3,
    "Edge Explorer": page4,
    "Trees": page5,
}


def show(page):
    return pages[page]()


starting_page = pn.state.session_args.get("page", [b"Overview"])[0].decode()
page = pn.widgets.RadioButtonGroup(
    value=starting_page,
    options=list(pages.keys()),
    name="Page",
    # sizing_mode="fixed",
    button_type="success",
    orientation="vertical",
)
ishow = pn.bind(show, page=page)
pn.state.location.sync(page, {"value": "page"})

ACCENT_COLOR = "#0072B5"
DEFAULT_PARAMS = {
    "site": "QC dashboard",
    "accent_base_color": ACCENT_COLOR,
    "header_background": ACCENT_COLOR,
}
pn.template.FastListTemplate(
    title=f"{trees_file}",
    sidebar=[page],
    main=[ishow],
    **DEFAULT_PARAMS,
).servable()
