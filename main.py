import pathlib
import time
import traceback

import click
import daiquiri
import tskit
import tszip

# Need to import daiquiri and set up logging before importing panel
# and bokeh, so we can set logging up correctly
daiquiri.setup(level="WARN")  # noqa
import panel as pn  # noqa

import model  # noqa
import pages  # noqa


logger = daiquiri.getLogger("tsqc")


def load_data(path):
    logger.info(f"Loading {path}")
    try:
        ts = tskit.load(path)
    except tskit.FileFormatError:
        ts = tszip.decompress(path)

    tsm = model.TSModel(ts, path.name)
    return tsm


page_map = {
    "Overview": pages.overview,
    "Mutations": pages.mutations,
    "Edges": pages.edges,
    "Edge Explorer": pages.edge_explorer,
    "Trees": pages.trees,
    "Nodes": pages.nodes,
    "Popgen": pages.popgen,
}


def get_app(tsm):
    pn.extension(sizing_mode="stretch_width")
    pn.extension("tabulator")

    def show(page):
        logger.info(f"Showing page {page}")
        yield pn.indicators.LoadingSpinner(value=True, width=50, height=50)
        try:
            before = time.time()
            content = page_map[page](tsm)
            duration = time.time() - before
            logger.info(f"Loaded page {page} in {duration:.2f}s")
        except Exception as e:
            error_message = f"An error occurred: {str(e)}"
            error_traceback = traceback.format_exc().replace("\n", "<br>")
            error_traceback = f"<pre>{error_traceback}</pre>"
            error_panel = pn.pane.Markdown(
                f"## Error\n\n{error_message}\n\n{error_traceback}",
                style={"color": "red"},
            )
            yield error_panel
            return
        yield content

    starting_page = pn.state.session_args.get("page", [b"Overview"])[0].decode()
    page = pn.widgets.RadioButtonGroup(
        value=starting_page,
        options=list(page_map.keys()),
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

    return pn.template.FastListTemplate(
        title=tsm.name,
        sidebar=[page],
        main=[ishow],
        **DEFAULT_PARAMS,
    )


def setup_logging(log_level, no_log_filter):
    if no_log_filter:
        logger = daiquiri.getLogger("root")
        logger.setLevel(log_level)
    else:
        # FIXME we can remove the "model" and "pages" bit when we've
        # structured as a package, as these should all have the prefix
        # tsqc
        # TODO figure out what's useful for users to track here, including
        # bokeh and tornado for now for dev
        loggers = ["tsqc", "model", "pages", "bokeh", "tornado"]
        for logname in loggers:
            logger = daiquiri.getLogger(logname)
            logger.setLevel(log_level)

        # Suppress traceback messages in the log like:
        # bokeh.core.serialization.DeserializationError: can't resolve reference 'p2091'
        # These pop up when the user clicks around between screens quickly,
        # and seems to occur when leaving a screen before it's fully rendered.
        # Some googling indicates that the messages are harmless, so simplest
        # to filter.
        logger = daiquiri.getLogger("bokeh.server.protocol_handler")
        logger.setLevel("CRITICAL")


@click.command()
@click.argument("path", type=click.Path(exists=True, dir_okay=False))
@click.option("--port", default=8080, help="Port to serve on")
@click.option(
    "--show/--no-show",
    default=True,
    help="Launch a web-browser showing the app",
)
@click.option("--log-level", default="INFO", help="Logging level")
@click.option(
    "--no-log-filter",
    default=False,
    is_flag=True,
    help="Do not filter the output log (advanced debugging only)",
)
def main(path, port, show, log_level, no_log_filter):
    """
    Run the tsqc server.
    """
    setup_logging(log_level, no_log_filter)
    tsm = load_data(pathlib.Path(path))

    # Note: functools.partial doesn't work here
    def app():
        return get_app(tsm)

    logger.info("Starting panel server")
    pn.serve(app, port=port, show=show, verbose=False)


if __name__ == "__main__":
    main()
