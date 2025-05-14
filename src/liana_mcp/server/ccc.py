from fastmcp import FastMCP, Context
import liana as li
import inspect
from pathlib import Path
import os
from ..schema.ccc import *
from scmcp_shared.util import add_op_log, savefig, filter_args, forward_request
from scmcp_shared.logging_config import setup_logger


ccc_mcp = FastMCP("LianaMCP-CCC-Server")

logger = setup_logger()

@ccc_mcp.tool()
async def ls_ccc_method():
    """List cell-cell communication method."""
    return str(li.mt.show_methods())


@ccc_mcp.tool()
async def communicate(
    request: CCCModel, 
    ctx: Context,
    sampleid: str = Field(default=None, description="adata sampleid for analysis"),
    dtype: str = Field(default="exp", description="the datatype of anndata.X"),
):
    """Cell-cell communication analysis with one method (cellphonedb, cellchat, connectome, natmi, etc.)"""

    try:
        result = await forward_request("ccc_communicate", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        method = request.method
        method_func = getattr(li.mt, method)
        func_kwargs = filter_args(request, method_func)
        method_func(adata, **func_kwargs)
        add_op_log(adata, method_func, func_kwargs)
        return [
            {"sampleid": sampleid or ads.active_id, "adata": adata, "dtype": dtype},
        ]
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e


@ccc_mcp.tool()
async def rank_aggregate(
    request: RankAggregateModel, 
    ctx: Context,
    sampleid: str = Field(default=None, description="adata sampleid for analysis"),
    dtype: str = Field(default="exp", description="the datatype of anndata.X"),
):
    """Get an aggregate of ligand-receptor scores from multiple Cell-cell communication methods."""

    try:
        result = await forward_request("ccc_rank_aggregate", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(dtype=dtype, sampleid=sampleid)
        func_kwargs = filter_args(request, li.mt.rank_aggregate)
        li.mt.rank_aggregate(adata, **func_kwargs)
        add_op_log(adata, li.mt.rank_aggregate, func_kwargs)
        return [
            {"sampleid": sampleid or ads.active_id, "adata": adata, "dtype": dtype},
        ]
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e
