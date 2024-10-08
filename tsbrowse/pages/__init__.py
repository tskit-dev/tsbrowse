from . import edges  # noqa: F401
from . import mutations  # noqa: F401
from . import nodes  # noqa: F401
from . import overview  # noqa: F401
from . import tables  # noqa: F401
from . import trees  # noqa: F401

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
