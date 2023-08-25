import holoviews as hv
import holoviews.operation.datashader as hd
import hvplot.pandas  # noqa
import numpy as np
import panel as pn

import config
from plot_helpers import filter_points
from plot_helpers import hover_points
from plot_helpers import make_hist


def page(tsm):
    hv.extension("bokeh")
    df_nodes = tsm.nodes_df
    df_internal_nodes = df_nodes[
        (df_nodes.is_sample == 0) & (df_nodes.ancestors_span != -np.inf)
    ]
    bins = min(50, int(np.sqrt(len(df_internal_nodes))))

    log_y_checkbox = pn.widgets.Checkbox(name="log y-axis of histogram", value=True)
    plot_options = pn.Column(
        pn.pane.Markdown("# Plot Options"),
        log_y_checkbox,
    )

    def make_node_hist_panel(tsm, log_y):
        nodes_hist = make_hist(
            df_internal_nodes.ancestors_span,
            "Ancestor spans per node",
            bins,
            log_y=log_y,
            plot_width=config.PLOT_WIDTH,
        )

        return pn.Column(pn.Row(nodes_hist))

    hist_panel = pn.bind(make_node_hist_panel, log_y=log_y_checkbox, tsm=tsm)

    points = df_nodes.hvplot.scatter(
        x="ancestors_span",
        y="time",
        hover_cols=["ancestors_span", "time"],
    ).opts(width=config.PLOT_WIDTH, height=config.PLOT_HEIGHT)

    range_stream = hv.streams.RangeXY(source=points)
    streams = [range_stream]
    filtered = points.apply(filter_points, streams=streams)
    hover = filtered.apply(hover_points, threshold=config.THRESHOLD)
    shaded = hd.datashade(filtered, width=400, height=400, streams=streams)

    main = (shaded * hover).opts(
        hv.opts.Points(tools=["hover"], alpha=0.1, hover_alpha=0.2, size=10)
    )

    plot_options = pn.Column(
        pn.pane.Markdown("# Plot Options"),
        log_y_checkbox,
    )

    anc_span_data = tsm.compute_ancestor_spans_heatmap_data()
    heatmap = hv.HeatMap(anc_span_data).opts(
        width=config.PLOT_WIDTH,
        height=config.PLOT_HEIGHT,
        tools=["hover"],
        colorbar=True,
    )

    return pn.Column(main, hist_panel, heatmap, plot_options)
