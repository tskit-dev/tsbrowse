import numpy as np
import panel as pn
import hvplot.pandas
import holoviews as hv
import pandas as pd
import holoviews.operation.datashader as hd
from holoviews.streams import RangeXY


# import hvplot.pandas
# from bokeh.sampledata.autompg import autompg
import tszip
import tskit
import utils


def mutations_data(ts):
    data = {
        "position": ts.sites_position[ts.mutations_site].astype(int),
        "node": ts.mutations_node,
        "time": ts.mutations_time,
    }
    return pd.DataFrame(data)


path = "/home/jk/work/github/sc2ts-paper/data/upgma-mds-1000-md-30-mm-3-2022-06-30-recinfo2-gisaid-il.ts.tsz"
ts = tszip.decompress(path)
# path = "/home/jk/work/github/sc2ts/results/full-md-30-mm-3-2021-04-13.ts"
# ts = tskit.load(path)
ti = utils.TreeInfo(ts, 1)

df_mutations = mutations_data(ts)


def page1():
    return pn.pane.HTML(ts)
    # hv_layout


def page2():
    points = df_mutations.hvplot.scatter("position", "time")
    pts = hd.datashade(points)  # , width=400, height=400)

    rasterized = hd.rasterize(points).opts(tools=["hover"])
    agg = hd.aggregate(points, width=120, height=120, streams=[RangeXY])
    dynamic = hv.util.Dynamic(agg, operation=hv.QuadMesh).opts(
        tools=["hover"], alpha=0, hover_alpha=0.2
    )
    spread = hd.dynspread(pts, threshold=0.8, how="over", max_px=5)
    main = (spread * dynamic).opts(width=1200)

    # FIXME need to link up the histogram here with RangeXY somehow
    ds_mutations = hv.Dataset(df_mutations)
    time_hist = hv.operation.element.histogram(ds_mutations, dimension="time")
    site_hist = hv.operation.element.histogram(ds_mutations, dimension="position")

    # stream = hv.streams.Tap(source=points, x=np.nan, y=np.nan)

    # @pn.depends(stream.param.x, stream.param.y)
    # def location(x, y):
    #     print("TAP", x, y)
    #     # This is a bad way to do it!
    #     # jitter = 10
    #     # pos = df_mutations.position.between(x - jitter, x + jitter)
    #     # time = df_mutations.time.between(y - jitter, y + jitter)
    #     # subset = df_mutations[pos & time]
    #     # print(subset)
    #     return pn.pane.Str(f"Click at {x:.2f}, {y:.2f}", width=200)

    count, bins = np.histogram(ti.sites_num_mutations, bins=range(29))
    overall_site_hist = hv.Histogram((count, bins)).opts(
        title="Mutations per site", tools=["hover"]
    )

    # Gah - these two axes are linked
    count, bins = np.histogram(ti.nodes_num_mutations, bins=range(10))
    overall_node_hist = hv.Histogram((count, bins)).opts(
        title="Mutations per node", tools=["hover"]
    )

    return pn.Column(
        main << time_hist << site_hist, overall_site_hist, overall_node_hist
    )


pn.extension(sizing_mode="stretch_width")

pages = {"Overview": page1, "Mutations": page2}


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
