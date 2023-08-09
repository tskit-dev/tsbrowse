import sys
import numpy as np
import panel as pn
import hvplot.pandas
import holoviews as hv
import pandas as pd
import holoviews.operation.datashader as hd
import tskit
import utils
import bokeh.models as bkm

# Usage: panel serve app.py --args /path/to/trees-file
path = sys.argv[1]
ts = tskit.load(path)
ti = utils.TreeInfo(ts, 1)


df_mutations = ti.mutations_data()
df_edges = ti.edges_data()


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


pn.extension(sizing_mode="stretch_width")
pn.extension("tabulator")

pages = {"Overview": page1, "Mutations": page2, "Edges": page3, "Edge Explorer": page4}


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
