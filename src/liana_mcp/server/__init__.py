from scmcp_shared.mcp_base import BaseMCPManager
from .preset.ccc import ccc_mcp
from .preset.pl import pl_mcp
from scmcp_shared.server.code import nb_mcp
from scmcp_shared.server.preset import io_mcp, ScanpyUtilMCP



ul_mcp = ScanpyUtilMCP(
    include_tools=["query_op_log", "check_samples"],
).mcp


class LianaMCPManager(BaseMCPManager):
    """Manager class for Liana MCP modules."""
    
    def init_mcp(self):
        """Initialize available Liana MCP modules."""
        self.available_modules = {
            "nb": nb_mcp,
            "ccc": ccc_mcp,
            "pl": pl_mcp,
            "io": io_mcp,
            "ul": ul_mcp,
        }
