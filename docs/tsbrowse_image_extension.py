import re
import time
from pathlib import Path

import panel as pn
from playwright.sync_api import sync_playwright
from sphinx.util import logging

from tsbrowse import app
from tsbrowse import model

logger = logging.getLogger(__name__)


class TSBrowsePreprocessor:
    def __init__(self, app):
        self.app = app
        self.src_dir = Path(app.srcdir)
        self.image_dir = self.src_dir / "_autoimages"
        self.image_dir.mkdir(parents=True, exist_ok=True)
        self.processed_images = set()

    def extract_image_specs(self, content):
        pattern = r"!\[(.*?)\]\(tsbrowse:(.*?):(.*?)\)"
        return re.findall(pattern, content)

    def generate_image(self, ts_file, page_name):
        output_filename = f"{Path(ts_file).stem}_{page_name}.png"
        output_path = self.image_dir / output_filename

        if output_path in self.processed_images:
            return output_filename

        logger.info(f"Generating screenshot for {ts_file}:{page_name}")

        tsm = model.TSModel(ts_file)
        app_ = app.App(tsm)
        url = "http://localhost:11337"
        server = pn.serve(app_.view, port=11337, threaded=True, show=False)
        try:
            # Wait for server to come up
            time.sleep(2)
            with sync_playwright() as p:
                browser = p.chromium.launch()
                page = browser.new_page()
                page.set_viewport_size({"width": 1920, "height": 1080})
                page.goto(url)
                page.get_by_role("button", name=page_name.title()).click()
                # Wait for page to load, would be better to have an element to wait for
                # But hard to do generically
                time.sleep(4)
                logger.info(f"Screenshotting to {output_path}")
                page.screenshot(path=str(output_path))
                browser.close()
        except Exception as e:
            logger.error(f"Failed to generate screenshot: {e}")
            return None
        finally:
            server.stop()

        self.processed_images.add(output_path)
        return output_filename

    def process_content(self, content, docname):
        image_specs = self.extract_image_specs(content)
        if not image_specs:
            return content

        modified_content = content
        for alt_text, ts_file, page_name in image_specs:
            image_filename = self.generate_image(ts_file, page_name)
            if image_filename:
                old_ref = f"![{alt_text}](tsbrowse:{ts_file}:{page_name})"
                new_ref = f"![{alt_text}](/_autoimages/{image_filename})"
                modified_content = modified_content.replace(old_ref, new_ref)

        return modified_content


def setup(app):
    logger.info("Setting up TSBrowse preprocessor")
    preprocessor = TSBrowsePreprocessor(app)

    def process_source(app, docname, source):
        content = source[0]
        if "tsbrowse:" not in content:
            return
        logger.info(f"Found tsbrowse references in {docname}")
        modified_content = preprocessor.process_content(content, docname)
        if modified_content != content:
            source[0] = modified_content

    app.connect("source-read", process_source)

    return {
        "version": "0.1",
        "parallel_read_safe": True,
        "parallel_write_safe": True,
    }
