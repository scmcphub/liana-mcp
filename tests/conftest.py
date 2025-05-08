
import pytest

@pytest.fixture
def mcp_config():
    return {
        "mcpServers": {
            "liana-mcp": {
                "command": "liana-mcp",
                "args": ["run"]
            }
        }
    }