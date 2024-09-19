import panel
import pytest
import tskit
import tszip

from tests import test_preprocess
from tsbrowse import model
from tsbrowse import pages
from tsbrowse import preprocess

# TODO give these some pytest metadata so they are named.
examples = [
    # No sites
    tskit.Tree.generate_balanced(5).tree_sequence,
    test_preprocess.single_tree_example_ts(),
    test_preprocess.single_tree_recurrent_mutation_example_ts(),
    test_preprocess.multiple_trees_example_ts(),
    test_preprocess.single_tree_with_polytomies_example_ts(),
]


class TestPages:
    @pytest.mark.parametrize("ts", examples)
    @pytest.mark.parametrize("page", pages.PAGES)
    def test_is_panel_layout_instance(self, ts, page, tmpdir):
        tszip.compress(ts, tmpdir / "test.tszip")
        preprocess.preprocess(tmpdir / "test.tszip", tmpdir / "test.tsbrowse")
        tsm = model.TSModel(tmpdir / "test.tsbrowse")
        page = page(tsm)
        assert isinstance(page.content, panel.layout.base.Panel)
        assert isinstance(page.sidebar, panel.layout.base.Panel)
