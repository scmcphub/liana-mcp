
import pytest

@pytest.fixture
def mcp():
    from liana_mcp.server import LianaMCPManager
    return LianaMCPManager("liana-mcp").mcp