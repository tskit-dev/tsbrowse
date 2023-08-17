import sys
import logging

import numpy as np
import panel as pn
import hvplot.pandas
import holoviews as hv
import pandas as pd

import tskit
import utils

import pathlib
import functools
import model
import pages


logger = logging.Logger(__file__)

# Usage: panel serve app.py --args /path/to/trees-file
path = pathlib.Path(sys.argv[1])
tsm = model.TSModel(tskit.load(path), path.name)


pn.extension(sizing_mode="stretch_width")
pn.extension("tabulator")

pages = {
    "Overview": pages.overview,
    "Mutations": pages.mutations,
    "Edges": pages.edges,
    "Edge Explorer": pages.edge_explorer,
    "Trees": pages.trees,
}


def show(page):
    return pages[page](tsm)


starting_page = pn.state.session_args.get("page", [b"Overview"])[0].decode()
page = pn.widgets.RadioButtonGroup(
    value=starting_page,
    options=list(pages.keys()),
    name="Page",
    # sizing_mode="fixed",
    button_type="success",
    orientation="vertical",
)
ishow = pn.bind(show, page=page)
pn.state.location.sync(page, {"value": "page"})

ACCENT_COLOR = "#0072B5"
DEFAULT_PARAMS = {
    "site": "QC dashboard",
    "accent_base_color": ACCENT_COLOR,
    "header_background": ACCENT_COLOR,
}
pn.template.FastListTemplate(
    title=tsm.name,
    sidebar=[page],
    main=[ishow],
    **DEFAULT_PARAMS,
).servable()
