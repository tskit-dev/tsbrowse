import pathlib
import sys
import time
import traceback

import daiquiri
import panel as pn
import tskit

import model
import pages


daiquiri.setup(level="INFO")
logger = daiquiri.getLogger("app")


@pn.cache()
def load_data(path):
    logger.info(f"Loading {path}")
    ts = tskit.load(path)
    tsm = model.TSModel(ts, path.name)
    return tsm


# Usage: panel serve app.py --args /path/to/trees-file
path = pathlib.Path(sys.argv[1])
tsm = load_data(path)

pn.extension(sizing_mode="stretch_width")
pn.extension("tabulator")

pages = {
    "Overview": pages.overview,
    "Mutations": pages.mutations,
    "Edges": pages.edges,
    "Edge Explorer": pages.edge_explorer,
    "Trees": pages.trees,
    "Nodes": pages.nodes,
    "Popgen": pages.popgen,
}


def show(page):
    logger.info(f"Showing page {page}")
    yield pn.indicators.LoadingSpinner(value=True, width=50, height=50)
    try:
        before = time.time()
        content = pages[page](tsm)
        duration = time.time() - before
        logger.info(f"Loaded page {page} in {duration:.2f}s")
    except Exception as e:
        error_message = f"An error occurred: {str(e)}"
        error_traceback = traceback.format_exc().replace("\n", "<br>")
        error_traceback = f"<pre>{error_traceback}</pre>"
        error_panel = pn.pane.Markdown(
            f"## Error\n\n{error_message}\n\n{error_traceback}", style={"color": "red"}
        )
        yield error_panel
        return
    yield content


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
