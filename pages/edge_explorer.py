import bokeh.models as bkm
import holoviews as hv
import panel as pn

import config


def page(tsm):
    hv.extension("bokeh")
    edges_df = tsm.edges_df
    node_id_input = pn.widgets.TextInput(value="", name="Node ID")
    edges_df["parent_time_right"] = edges_df["parent_time"]
    tabulator = pn.widgets.Tabulator(show_index=False)

    def plot_data(node_id):
        if len(node_id) > 0:
            filtered_df = edges_df[edges_df["child"] == int(node_id)]
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
                width=config.PLOT_WIDTH,
                height=config.PLOT_HEIGHT,
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
