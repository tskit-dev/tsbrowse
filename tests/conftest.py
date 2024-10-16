import panel as pn
import pytest

PORT = [6000]


@pytest.fixture
def port():
    PORT[0] += 1
    return PORT[0]


@pytest.fixture(autouse=True)
def server_cleanup():
    """
    Clean up server state after each test.
    """
    try:
        yield
    finally:
        pn.state.reset()


def pytest_addoption(parser):
    parser.addoption(
        "--save-screenshots",
        action="store_true",
        default=False,
        help="Save screenshots during tests",
    )
