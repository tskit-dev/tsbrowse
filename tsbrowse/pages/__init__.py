from . import (
    edges,  # noqa: F401
    mutations,  # noqa: F401
    nodes,  # noqa: F401
    overview,  # noqa: F401
    tables,  # noqa: F401
    trees,  # noqa: F401
)

PAGES = [
    overview.OverviewPage,
    tables.TablesPage,
    mutations.MutationsPage,
    edges.EdgesPage,
    trees.TreesPage,
    nodes.NodesPage,
]
PAGES_MAP = {page.key: page for page in PAGES}
PAGES_BY_TITLE = {page.title: page for page in PAGES}
