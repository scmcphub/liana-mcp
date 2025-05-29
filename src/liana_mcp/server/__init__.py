from scmcp_shared.server import BaseMCPManager
from scmcp_shared.server import io_mcp, ScanpyUtilMCP
from .ccc import ccc_mcp
from .pl import pl_mcp


ul_mcp = ScanpyUtilMCP(
    include_tools=["query_op_log", "check_samples"],
).mcp



class LianaMCPManager(BaseMCPManager):
    """Manager class for Liana MCP modules."""
    
    def _init_modules(self):
        """Initialize available Liana MCP modules."""
        self.available_modules = {
            "io": io_mcp,
            "ccc": ccc_mcp,
            "pl": pl_mcp,
            "ul": ul_mcp
        }
