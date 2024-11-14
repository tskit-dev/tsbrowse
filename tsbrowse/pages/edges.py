import bokeh.models as bkm
import holoviews as hv
import holoviews.operation.datashader as hd
import numpy as np
import panel as pn

from .. import config
from ..plot_helpers import filter_points
from ..plot_helpers import hover_points
from ..plot_helpers import make_hist


def parse_range(range_str):
    try:
        start, end = map(float, range_str.split(":"))
        return (start, end)
    except (ValueError, TypeError):
        return None  # Default range if parsing fails


def make_edges_panel(log_y, node_type, x_range, tsm):
    edges_df = tsm.edges_df
    if node_type == "Child node":
        time_column = "child_time"
        y_label = f"Time of Child Node ({tsm.ts.time_units})"
    else:
        time_column = "parent_time"
        y_label = f"Time of Parent Node ({tsm.ts.time_units})"

    if log_y:
        log_time_column = f"log_{time_column}"
        edges_df[log_time_column] = np.log10(1 + edges_df[time_column])
        y_dim_left = log_time_column
        y_dim_right = log_time_column + "_right"
        y_label = f"log ( {y_label} )"
    else:
        y_dim_left = time_column
        y_dim_right = time_column + "_right"
        y_label = y_label

    edges_df[y_dim_right] = edges_df[y_dim_left]

    lines = hv.Segments(edges_df, kdims=["left", y_dim_left, "right", y_dim_right])
    range_stream = hv.streams.RangeXY(source=lines)
    streams = [range_stream]
    filtered = lines.apply(filter_points, streams=streams)
    hover = filtered.apply(hover_points)
    shaded = hd.datashade(lines, streams=streams, cmap=config.PLOT_COLOURS[1:])
    hover_tool = bkm.HoverTool(
        tooltips=[
            ("child", "@child"),
            ("parent", "@parent"),
            ("span", "@span"),
            ("branch_length", "@branch_length"),
        ]
    )
    opts = {
        "tools": [hover_tool],
        "width": config.PLOT_WIDTH,
        "height": config.PLOT_HEIGHT,
        "xlabel": "Position",
        "ylabel": y_label,
    }
    x_range = parse_range(x_range)
    if x_range is not None:
        opts["xlim"] = x_range
    main = (shaded * hover).opts(
        hv.opts.Segments(
            **opts,
        )
    )

    edges_df["time_span"] = edges_df["parent_time"] - edges_df["child_time"]
    gspan_hist = make_hist(
        edges_df["span"],
        "Genomic span",
        "auto",
        log_y=log_y,
        xlabel="Genomic Span",
        ylabel="Number of Edges",
    )
    tspan_hist = make_hist(
        edges_df["time_span"],
        "Time span",
        "auto",
        log_y=log_y,
        xlabel=f"Time span ({tsm.ts.time_units})",
        ylabel="Number of Edges",
    )
    area_hist = make_hist(
        edges_df["span"] * edges_df["time_span"],
        "Edge area",
        "auto",
        log_y=log_y,
        xlabel="Time Span * Genomic Span",
        ylabel="Number of Edges",
    )

    area_hist.opts(width=None, height=None, responsive=True, xrotation=90)
    gspan_hist.opts(width=None, height=None, responsive=True, xrotation=90)
    tspan_hist.opts(width=None, height=None, responsive=True, xrotation=90)

    main = main.opts(width=None, height=None, responsive=True)

    gspec = pn.GridSpec(sizing_mode="stretch_both", ncols=5, nrows=3)

    gspec[:, :3] = main
    gspec[0, 3:] = gspan_hist
    gspec[1, 3:] = tspan_hist
    gspec[2, 3:] = area_hist

    return gspec


class EdgesPage:
    key = "edges"
    title = "Edges"

    def __init__(self, tsm):
        log_y_checkbox = pn.widgets.Checkbox(name="Log Y-axis", value=False)
        node_type_radio = pn.widgets.RadioBoxGroup(
            options=["Parent node", "Child node"], value="Parent node", inline=True
        )
        # using a markdown widget to display radiobox title
        # (https://github.com/holoviz/panel/issues/1313):
        radio_title = pn.pane.Markdown("Plot time of:")
        x_range = pn.widgets.TextInput(value="", name="X Range (start:stop)")
        options_box = pn.WidgetBox(
            log_y_checkbox, radio_title, node_type_radio, x_range
        )
        edges_panel = pn.bind(
            make_edges_panel,
            log_y=log_y_checkbox,
            node_type=node_type_radio,
            x_range=x_range,
            tsm=tsm,
        )
        self.content = pn.Column(
            edges_panel,
            sizing_mode="stretch_both",
        )
        self.sidebar = pn.Column(pn.pane.Markdown("# Edges"), options_box)
