from pathlib import Path

import click
import daiquiri

# Need to import daiquiri and set up logging before importing panel
# and bokeh, so we can set logging up correctly
daiquiri.setup(level="WARN")  # noqa
import panel as pn  # noqa

from . import model  # noqa
from . import pages  # noqa
from . import config  # noqa
from . import preprocess as preprocess_  # noqa
from . import raster  # noqa
from . import app  # noqa


logger = daiquiri.getLogger("tsbrowse")


def setup_logging(log_level, no_log_filter):
    if no_log_filter:
        logger = daiquiri.getLogger("root")
        logger.setLevel(log_level)
    else:
        # TODO figure out what's useful for users to track here, including
        # bokeh and tornado for now for dev
        loggers = ["tsbrowse", "cache", "bokeh", "tornado"]
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


@click.group()
def cli():
    """Command line interface for tsbrowse."""
    pass


@cli.command()
@click.argument("path", type=click.Path(exists=True, dir_okay=False))
@click.option("--annotations-file", type=click.Path(exists=True, dir_okay=False))
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
def serve(path, port, show, log_level, no_log_filter, annotations_file):
    """
    Run the tsbrowse server.
    """
    setup_logging(log_level, no_log_filter)

    tsm = model.TSModel(path)
    if annotations_file:
        config.ANNOTATIONS_FILE = annotations_file

    logger.info("Starting panel server")
    app_ = app.App(tsm)
    pn.serve(app_.view, port=port, show=show, verbose=False)


@cli.command()
@click.argument("tszip_path", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--output",
    "-o",
    type=click.Path(dir_okay=False),
    default=None,
    help="Optional output filename, defaults to tszip_path with .tsbrowse extension",
)
@click.option(
    "-v",
    "--verbose",
    count=True,
    help="Increase verbosity. Use -v for INFO and -vv for DEBUG.",
)
def preprocess(tszip_path, output, verbose):
    """
    Preprocess a tskit tree sequence or tszip file, producing a .tsbrowse file.
    """
    tszip_path = Path(tszip_path)
    if output is None:
        output = tszip_path.with_suffix(".tsbrowse")

    # Set log level based on verbosity
    if verbose == 1:
        log_level = "INFO"
    elif verbose >= 2:
        log_level = "DEBUG"
    else:
        log_level = "WARN"

    setup_logging(log_level, no_log_filter=False)

    preprocess_.preprocess(tszip_path, output, show_progress=True)
    logger.info(f"Preprocessing completed. Output saved to: {output}")
    print(f"Preprocessing completed. You can now view with `tsbrowse serve {output}`")


@cli.command()
@click.argument("tsbrowse_path", type=click.Path(exists=True, dir_okay=False))
@click.argument(
    "page",
    type=click.Choice([page.key for page in pages.PAGES]),
)
@click.option(
    "--output",
    "-o",
    type=click.Path(dir_okay=False),
    default=None,
    help="Optional output filename for the screenshot. If not provided, it will"
    "be automatically generated based on the input filename and page name.",
)
def screenshot(tsbrowse_path, page, output):
    """
    Create a screenshot of a specific page from a .tsbrowse file.

    TSBROWSE_PATH: Path to the .tsbrowse file to screenshot.

    PAGE: The specific page to capture (e.g., 'overview', 'mutations', 'edges').

    Example usage:
    tsbrowse screenshot example.tsbrowse overview --output outfile.png
    """
    if output is None:
        base_name = tsbrowse_path.rsplit(".", 1)[0]
        output = f"{base_name}_{page}.png"
    tsm = model.TSModel(tsbrowse_path)
    page = pages.PAGES_MAP[page](tsm)
    raster.raster_component(page, output)
    logger.info(f"Screenshot saved to: {output}")


if __name__ == "__main__":
    cli()
