import time

import msprime
import panel as pn
import pytest
from playwright.sync_api import expect

from tsbrowse import app
from tsbrowse import model
from tsbrowse import preprocess


@pytest.fixture
def save_screenshots(request):
    return request.config.getoption("--save-screenshots")


def test_component(page, port, tmpdir, save_screenshots):
    ts = msprime.sim_ancestry(
        10, sequence_length=1000, recombination_rate=1e-2, random_seed=1
    )
    ts = msprime.sim_mutations(ts, rate=1e-2, random_seed=1)
    ts.dump(tmpdir / "test.trees")
    preprocess.preprocess(tmpdir / "test.trees", tmpdir / "test.tsbrowse")

    tsm = model.TSModel(tmpdir / "test.tsbrowse")
    component = app.App(tsm)

    url = f"http://localhost:{port}"
    server = pn.serve(component.view, port=port, threaded=True, show=False)
    time.sleep(2)
    page.goto(url)

    page.set_viewport_size({"width": 1920, "height": 1080})
    expect(page.get_by_role("link", name="Tree Sequence")).to_be_visible()
    expect(page.get_by_role("cell", name="Provenance Timestamp")).to_be_visible()
    if save_screenshots:
        page.screenshot(path="overview.png")

    page.get_by_role("button", name="Tables").click()
    expect(page.get_by_label("Select Table")).to_be_visible()
    expect(page.get_by_placeholder("Enter query expression (e.g")).to_be_visible()
    page.get_by_label("Select Table").select_option("trees")
    expect(page.get_by_text("total_branch_length")).to_be_visible()
    if save_screenshots:
        page.screenshot(path="tables.png")

    page.get_by_role("button", name="Mutations").click()
    expect(page.get_by_text("Log Y-axis")).to_be_visible()
    expect(page.get_by_title("Reset").locator("div")).to_be_visible()
    if save_screenshots:
        page.screenshot(path="mutations.png")
    # This test was flakey on CI, despite the fact that it works locally.
    # This horrendous selector is needed because the click event is not
    # handled be the canvas, but by a floating div on top.
    # page.locator("div:nth-child(6) > .bk-Canvas > div:nth-child(12)").click(
    #     position={"x": 558, "y": 478}
    # )
    # expect(page.get_by_text("Mutation information")).to_be_visible()
    # expect(page.get_by_text("0.52")).to_be_visible()
    # if save_screenshots:
    #     page.screenshot(path="mutations-popup.png")

    page.get_by_role("button", name="Edges").click()
    expect(page.get_by_text("Parent node")).to_be_visible()
    expect(page.get_by_title("Reset").locator("div")).to_be_visible()
    if save_screenshots:
        page.screenshot(path="edges.png")

    page.get_by_role("button", name="Trees").click()
    expect(page.get_by_text("log y-axis")).to_be_visible()
    expect(page.locator(".bk-Canvas > div:nth-child(12)").first).to_be_visible()
    if save_screenshots:
        page.screenshot(path="trees.png")

    page.get_by_role("button", name="Nodes").click()
    expect(page.get_by_role("heading", name="Node Flags")).to_be_visible()
    expect(page.get_by_title("Reset").locator("div")).to_be_visible()
    if save_screenshots:
        page.screenshot(path="nodes.png")

    server.stop()
