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

        used_bits = np.uint32(0)
        flag_counts = np.zeros(32, dtype=int)
        for flag in df_nodes["flags"]:
            used_bits |= np.uint32(flag)
            for i in range(32):
                if flag & (1 << i):
                    flag_counts[i] += 1

        used_flags = [i for i in range(32) if used_bits & (1 << i)]
        options = {f"Flag {i} ({flag_counts[i]})": i for i in used_flags}
        node_flag_checkboxes = pn.widgets.CheckBoxGroup(
            name="Node Types", value=[], options=options
        )

        def make_node_hist_panel(log_y, node_flags):
            df = df_nodes[
                df_nodes["flags"].apply(
                    lambda x: all((x & (1 << flag)) != 0 for flag in node_flags)
                )
            ]
            nodes_hist = make_hist(
                df.ancestors_span,
                "Ancestor spans per node",
                bins,
                log_y=log_y,
                xlabel="Ancestor Span",
                ylabel="Number of Nodes",
            )
            nodes_hist.opts(width=None, height=None, responsive=True)
            return nodes_hist

        def make_node_plot(node_flags):
            df = df_nodes[
                df_nodes["flags"].apply(
                    lambda x: all((x & (1 << flag)) != 0 for flag in node_flags)
                )
            ]
            points = df.hvplot.scatter(
                x="ancestors_span",
                y="time",
                hover_cols=["id", "ancestors_span", "time"],
            ).opts(
                xlabel="Ancestor Span",
                ylabel=f"Time ({tsm.ts.time_units})",
            )

            range_stream = hv.streams.RangeXY(source=points)
            streams = [range_stream]
            filtered = points.apply(filter_points, streams=streams)
            hover = filtered.apply(hover_points, threshold=config.THRESHOLD)
            shaded = hd.datashade(
                points,
                streams=streams,
                cmap=config.PLOT_COLOURS[1:],
            )

            main = (shaded * hover).opts(
                hv.opts.Points(tools=["hover"], alpha=0.1, hover_alpha=0.2, size=10)
            )
            main.opts(width=None, height=None, responsive=True)
            return main

        gspec = pn.GridSpec(sizing_mode="stretch_both", ncols=3)
        gspec[0:2, 0] = pn.bind(make_node_plot, node_flags=node_flag_checkboxes)
        gspec[2:, 0] = pn.bind(
            make_node_hist_panel, log_y=log_y_checkbox, node_flags=node_flag_checkboxes
        )

        self.content = gspec
        self.sidebar = pn.Column(
            pn.pane.Markdown("# Nodes"),
            log_y_checkbox,
            pn.pane.Markdown(f"### Total Nodes: **{len(df_nodes)}**"),
            pn.pane.Markdown("### Node Flags"),
            node_flag_checkboxes,
        )
