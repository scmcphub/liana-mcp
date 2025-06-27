
import pytest

@pytest.fixture
def mcp():
    from liana_mcp.server import LianaMCPManager
    from scmcp_shared.backend import AdataManager
    return LianaMCPManager("liana-mcp", backend=AdataManager).mcp