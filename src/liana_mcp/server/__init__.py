import asyncio
from fastmcp import FastMCP
from collections.abc import AsyncIterator
from contextlib import asynccontextmanager
from typing import Any

import scmcp_shared.server as shs
from scmcp_shared.util import filter_tools

from .ccc import ccc_mcp
from.pl import pl_mcp


ads = shs.AdataState()

@asynccontextmanager
async def adata_lifespan(server: FastMCP) -> AsyncIterator[Any]:
    yield ads

liana_mcp = FastMCP("Liana-MCP-Server", lifespan=adata_lifespan)


async def setup(modules=None):
    ul_mcp = await filter_tools(shs.ul_mcp, include_tools=["query_op_log", "check_samples"])
    mcp_dic = {
        "io": shs.io_mcp, 
        "ccc": ccc_mcp, 
        "pl": pl_mcp, 
        "ul": ul_mcp
        }
    if modules is None or modules == "all":
        modules = ["io", "ccc", "pl", "ul"]
    for module in modules:
        await liana_mcp.import_server(module, mcp_dic[module])
