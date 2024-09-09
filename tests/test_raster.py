import msprime
import pytest

from tsbrowse import pages
from tsbrowse import raster


@pytest.fixture(scope="module")
def ts():
    ts = msprime.sim_ancestry(
        recombination_rate=1e-3,
        samples=100,
        sequence_length=10_000,
    )
    return msprime.sim_mutations(ts, rate=1e-2, random_seed=43)


class TestRaster:
    @pytest.mark.parametrize("page", pages.PAGES_MAP.values())
    def test_mutation_scatter(self, page, ts, tmp_path):
        print(ts)
        raster.raster_component(
            page.page,
            ts,
            tmp_path / "image.png",
        )
        assert (tmp_path / "image.png").exists()
