import bokeh.models as bkm
import holoviews as hv
import holoviews.operation.datashader as hd
import numpy as np
import panel as pn

from .. import config
from ..plot_helpers import filter_points
from ..plot_helpers import hover_points


def make_edges_panel(log_y, tsm):
    edges_df = tsm.edges_df
    if log_y:
        edges_df["log_parent_time"] = np.log10(1 + edges_df["parent_time"])
        edges_df["log_parent_time_right"] = edges_df["log_parent_time"]
        y_dim_left = "log_parent_time"
        y_dim_right = "log_parent_time_right"
    else:
        edges_df["parent_time_right"] = edges_df["parent_time"]
        y_dim_left = "parent_time"
        y_dim_right = "parent_time_right"

    lines = hv.Segments(edges_df, kdims=["left", y_dim_left, "right", y_dim_right])
    range_stream = hv.streams.RangeXY(source=lines)
    streams = [range_stream]
    filtered = lines.apply(filter_points, streams=streams)
    hover = filtered.apply(hover_points)
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
            width=config.PLOT_WIDTH,
            height=config.PLOT_HEIGHT,
            xlabel="Position",
            ylabel="Time (parent node)",
        )
    )
    return pn.Column(main)


def page(tsm):
    hv.extension("bokeh")

    log_y_checkbox = pn.widgets.Checkbox(name="Log y-axis", value=False)
    plot_options = pn.Column(
        pn.pane.Markdown("### Plot Options"),
        log_y_checkbox,
    )
    edges_panel = pn.bind(make_edges_panel, log_y=log_y_checkbox, tsm=tsm)
    return pn.Column(plot_options, edges_panel)
