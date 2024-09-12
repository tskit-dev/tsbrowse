from . import edge_explorer  # noqa: F401
from . import edges  # noqa: F401
from . import frequency_spectra  # noqa: F401
from . import mutations  # noqa: F401
from . import nodes  # noqa: F401
from . import overview  # noqa: F401
from . import popgen  # noqa: F401
from . import trees  # noqa: F401

PAGES_MAP = {
    "Overview": overview,
    "Mutations": mutations,
    "Edges": edges,
    "Edge Explorer": edge_explorer,
    "Trees": trees,
    "Nodes": nodes,
    "Popgen": popgen,
    "SFS": frequency_spectra,
}
