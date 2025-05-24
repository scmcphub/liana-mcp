import asyncio
import pytest
from fastmcp import Client
import anndata
from pathlib import Path
from liana_mcp.server import liana_mcp, setup


asyncio.run(setup())


@pytest.mark.asyncio 
async def test_ls_ccc_method():
    """Test listing available CCC methods."""
    async with Client(liana_mcp) as client:
        result = await client.call_tool("ccc_ls_ccc_method", {})
        assert isinstance(result[0].text, str)
        assert "cellphonedb" in result[0].text.lower()
        assert "cellchat" in result[0].text.lower()

@pytest.mark.asyncio 
async def test_ccc_communicate():
    """Test cell-cell communication analysis with different methods."""
    test_dir = Path(__file__).parent / "data/pbmc68k_reduced.h5ad"
    
    async with Client(liana_mcp) as client:
        # First read the data
        result = await client.call_tool("io_read", {"request": {"filename": test_dir}})
        assert "AnnData" in result[0].text

        # Test cellphonedb method
        result = await client.call_tool("ccc_communicate", {
            "request": {
                "method": "cellphonedb",
                "groupby": "bulk_labels"
            }
        })
        assert "adata" in result[0].text

        # Test cellchat method
        result = await client.call_tool("ccc_communicate", {
            "request": {
                "method": "cellchat",
                "groupby": "bulk_labels"
            }
        })
        assert "adata" in result[0].text

@pytest.mark.asyncio 
async def test_rank_aggregate():
    """Test rank aggregation of multiple CCC methods."""
    test_dir = Path(__file__).parent / "data/pbmc68k_reduced.h5ad"
    
    async with Client(liana_mcp) as client:
        # First read the data
        result = await client.call_tool("io_read", {"request": {"filename": test_dir}, "adinfo":{ "sampleid": "pbmc68k", "adtype": "exp"}})
        assert "AnnData" in result[0].text

        # Run rank aggregation
        result = await client.call_tool("ccc_rank_aggregate", {
            "request": {
                "methods": ["cellphonedb", "cellchat"],
                "groupby": "bulk_labels"
            }
        })
        assert "adata" in result[0].text
