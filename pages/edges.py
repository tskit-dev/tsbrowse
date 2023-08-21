import bokeh.models as bkm
import holoviews as hv
import holoviews.operation.datashader as hd
import panel as pn

import config
from plot_helpers import filter_points
from plot_helpers import hover_points


def page(tsm):
    hv.extension("bokeh")
    edges_df = tsm.edges_df
    edges_df["parent_time_right"] = edges_df["parent_time"]
    lines = hv.Segments(
        edges_df, kdims=["left", "parent_time", "right", "parent_time_right"]
    )
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
            ylabel="Time",
        )
    )

    return pn.Column(main)
