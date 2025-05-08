import pytest
from fastmcp import Client
from pathlib import Path

@pytest.mark.asyncio 
async def test_read_and_write(mcp_config):
    # Pass the server directly to the Client constructor
    test_dir = Path(__file__).parent / "data/pbmc68k_reduced.h5ad"
    outfile = Path(__file__).parent / "data/test.h5ad"
    async with Client(mcp_config) as client:
        result = await client.call_tool("io_read", {"request":{"filename": test_dir}})
        assert "AnnData" in result[0].text

        result = await client.call_tool("ccc_communicate", {"request":{"method": "cellphonedb", "groupby": "bulk_labels"}})
        result = await client.call_tool("pl_circle_plot", {"request":{"uns_key": "liana_res", "groupby": "bulk_labels"}})
        assert "figpath" in result[0].text

