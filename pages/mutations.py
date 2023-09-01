import holoviews as hv
import holoviews.operation.datashader as hd
import hvplot.pandas  # noqa
import panel as pn

import config
from plot_helpers import filter_points
from plot_helpers import hover_points
from plot_helpers import make_hist
from plot_helpers import make_hist_on_axis


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


def page(tsm):
    hv.extension("bokeh")
    plot_width = 1000
    log_y_checkbox = pn.widgets.Checkbox(
        name="Log y-axis of Mutations per site/node plots", value=False
    )

    points = tsm.mutations_df.hvplot.scatter(
        x="position",
        y="time",
        hover_cols=["position", "time", "mutation_node", "node_flag"],
    ).opts(width=plot_width, height=config.PLOT_HEIGHT)

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
    hover = filtered.apply(hover_points)
    shaded = hd.datashade(filtered, width=400, height=400, streams=streams)

    main = (shaded * hover).opts(
        hv.opts.Points(tools=["hover"], alpha=0.1, hover_alpha=0.2, size=10)
    )

    hist_panel = pn.bind(make_hist_panel, log_y=log_y_checkbox, tsm=tsm)

    plot_options = pn.Column(
        pn.pane.Markdown("## Plot Options"),
        log_y_checkbox,
    )

    return pn.Column(main << time_hist << site_hist, hist_panel, plot_options)
