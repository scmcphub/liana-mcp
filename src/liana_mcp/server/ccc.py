from fastmcp import FastMCP, Context
from fastmcp.exceptions import ToolError
import liana as li
import inspect
from pathlib import Path
import os
from ..schema.ccc import *
from scmcp_shared.util import add_op_log, savefig, filter_args, forward_request, get_ads
from scmcp_shared.logging_config import setup_logger
from scmcp_shared.schema import AdataModel

ccc_mcp = FastMCP("LianaMCP-CCC-Server")

logger = setup_logger()

@ccc_mcp.tool()
async def ls_ccc_method():
    """List cell-cell communication method."""
    return str(li.mt.show_methods())


@ccc_mcp.tool()
async def communicate(
    request: CCCModel,
    adinfo: AdataModel = AdataModel()
):
    """Cell-cell communication analysis with one method (cellphonedb, cellchat, connectome, natmi, etc.)"""

    try:
        result = await forward_request("ccc_communicate", request, adinfo)
        if result is not None:
            return result
        ads = get_ads()
        adata = ads.get_adata(adinfo=adinfo)
        method = request.method
        method_func = getattr(li.mt, method)
        func_kwargs = filter_args(request, method_func)
        method_func(adata, **func_kwargs)
        add_op_log(adata, method_func, func_kwargs, adinfo)
        return [
            {"sampleid": adinfo.sampleid or ads.active_id, "adtype": adinfo.adtype, "adata": adata},
        ]
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)



@ccc_mcp.tool()
async def rank_aggregate(
    request: RankAggregateModel, 
    adinfo: AdataModel = AdataModel()
):
    """Get an aggregate of ligand-receptor scores from multiple Cell-cell communication methods."""

    try:
        result = await forward_request("ccc_rank_aggregate", request, adinfo)
        if result is not None:
            return result
        ads = get_ads()
        adata = ads.get_adata(adinfo=adinfo)
        func_kwargs = filter_args(request, li.mt.rank_aggregate)
        li.mt.rank_aggregate(adata, **func_kwargs)
        add_op_log(adata, li.mt.rank_aggregate, func_kwargs, adinfo)
        return [
            {"sampleid": adinfo.sampleid or ads.active_id, "adtype": adinfo.adtype, "adata": adata},
        ]
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)
