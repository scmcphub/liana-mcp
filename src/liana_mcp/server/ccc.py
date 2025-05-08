from fastmcp import FastMCP, Context
import liana as li
import inspect
from pathlib import Path
import os
from ..schema.ccc import *
from ..util import add_op_log, savefig, filter_args, forward_request
from ..logging_config import setup_logger

ccc_mcp = FastMCP("LianaMCP-CCC-Server")

logger = setup_logger()

@ccc_mcp.tool()
async def ls_ccc_method(request: ListCCCMethodModel, ctx: Context):
    """List cell-cell communication method."""
    return str(li.mt.show_methods())


@ccc_mcp.tool()
async def communicate(request: CCCModel, ctx: Context):
    """Cell-cell communication analysis with one method (cellphonedb, cellchat, connectome, natmi, etc.)"""
    result = await forward_request("ccc_communicate", request.model_dump())
    if result is not None:
        return result
    ads = ctx.request_context.lifespan_context
    adata = ads.adata_dic[ads.active_id]
    method = request.method
    method_func = getattr(li.mt, method)
    func_kwargs = filter_args(request, method_func)
    method_func(adata, **func_kwargs)
    add_op_log(adata, method_func, func_kwargs)
    return adata


@ccc_mcp.tool()
async def rank_aggregate(request: RankAggregateModel, ctx: Context):
    """Get an aggregate of ligand-receptor scores from multiple Cell-cell communication methods."""
    result = await forward_request("ccc_rank_aggregate", request.model_dump())
    if result is not None:
        return result
    ads = ctx.request_context.lifespan_context
    adata = ads.adata_dic[ads.active_id]
    func_kwargs = filter_args(request, li.mt.rank_aggregate)
    li.mt.rank_aggregate(adata, **func_kwargs)
    add_op_log(adata, li.mt.rank_aggregate, func_kwargs)
    return adata
