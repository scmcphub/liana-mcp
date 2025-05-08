
import asyncio
from fastmcp import FastMCP
from .server import ccc_mcp, pl_mcp, io_mcp


mcp = FastMCP("Liana-MCP-Server")


async def setup():
    await mcp.import_server("io", io_mcp)
    await mcp.import_server("ccc", ccc_mcp)
    await mcp.import_server("pl", pl_mcp)



if __name__ == "__main__":
    asyncio.run(setup())
    mcp.run()
    