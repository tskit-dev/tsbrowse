import panel
import pytest
import tskit

from tests import test_data_model
from tsbrowse import model
from tsbrowse import pages

# TODO give these some pytest metadata so they are named.
examples = [
    # No sites
    tskit.Tree.generate_balanced(5).tree_sequence,
    test_data_model.single_tree_example_ts(),
    test_data_model.single_tree_recurrent_mutation_example_ts(),
    test_data_model.multiple_trees_example_ts(),
    test_data_model.single_tree_with_polytomies_example_ts(),
]

display_pages = [
    pages.overview,
    pages.mutations,
    pages.edges,
    pages.edge_explorer,
    pages.trees,
    pages.nodes,
    pages.popgen,
]


class TestPages:
    @pytest.mark.parametrize("ts", examples)
    @pytest.mark.parametrize("page", display_pages)
    def test_is_panel_layout_instance(self, ts, page):
        tsm = model.TSModel(ts)
        ui = page.page(tsm)
        assert isinstance(ui, panel.layout.base.Panel)
