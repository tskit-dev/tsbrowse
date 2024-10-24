import holoviews as hv
import holoviews.operation.datashader as hd
import hvplot.pandas  # noqa: F401
import numpy as np
import panel as pn

from .. import config
from ..plot_helpers import filter_points
from ..plot_helpers import hover_points
from ..plot_helpers import make_hist


class NodesPage:
    key = "nodes"
    title = "Nodes"

    def __init__(self, tsm):
        df_nodes = tsm.nodes_df
        df_nodes = df_nodes[(df_nodes.ancestors_span != -np.inf)]
        bins = min(50, int(np.sqrt(len(df_nodes))))

        log_y_checkbox = pn.widgets.Checkbox(name="Log Y-axis of Histogram", value=True)

        def make_node_hist_panel(tsm, log_y):
            nodes_hist = make_hist(
                df_nodes.ancestors_span,
                "Ancestor spans per node",
                bins,
                log_y=log_y,
                plot_width=config.PLOT_WIDTH,
                xlabel="Ancestor Span",
                ylabel="Number of Nodes",
            )

            return pn.Column(pn.Row(nodes_hist))

        hist_panel = pn.bind(make_node_hist_panel, log_y=log_y_checkbox, tsm=tsm)

        def make_node_plot(data, node_types):
            df = data[data["flags"].isin(node_types)]
            points = df.hvplot.scatter(
                x="ancestors_span",
                y="time",
                hover_cols=["id", "ancestors_span", "time"],
            ).opts(
                width=config.PLOT_WIDTH,
                height=config.PLOT_HEIGHT,
                xlabel="Ancestor Span",
                ylabel=f"Time ({tsm.ts.time_units})",
            )

            range_stream = hv.streams.RangeXY(source=points)
            streams = [range_stream]
            filtered = points.apply(filter_points, streams=streams)
            hover = filtered.apply(hover_points, threshold=config.THRESHOLD)
            shaded = hd.datashade(
                points,
                width=400,
                height=400,
                streams=streams,
                cmap=config.PLOT_COLOURS[1:],
            )

            main = (shaded * hover).opts(
                hv.opts.Points(tools=["hover"], alpha=0.1, hover_alpha=0.2, size=10)
            )
            return main

        def make_node_panel(node_types):
            nodes_spans_plot = make_node_plot(
                df_nodes,
                node_types=node_types,
            )
            return pn.Row(nodes_spans_plot)

        anc_options = list(df_nodes["flags"].unique())
        checkboxes = pn.widgets.CheckBoxGroup(
            name="Node Types", value=anc_options, options=anc_options
        )
        nodes_panel = pn.bind(make_node_panel, node_types=checkboxes)
        node_options = pn.Column(
            pn.pane.Markdown("### Node Flags"),
            checkboxes,
        )

        self.content = pn.Column(nodes_panel, hist_panel)
        self.sidebar = pn.Column(
            pn.pane.Markdown("# Nodes"), node_options, log_y_checkbox
        )
