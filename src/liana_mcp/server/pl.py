from fastmcp import FastMCP, Context
import liana as li
import inspect
from pathlib import Path
import os
from ..schema.pl import *
from ..util import add_op_log, savefig, filter_args, forward_request
from ..logging_config import setup_logger
logger = setup_logger()


pl_mcp = FastMCP("LianaMCP-CCC-Server")


@pl_mcp.tool()
async def circle_plot(request: CirclePlotModel, ctx: Context):
    """Visualize cell-cell communication network using a circular plot."""
    result = await forward_request("pl_circle_plot", request.model_dump())
    if result is not None:
        return result
    ads = ctx.request_context.lifespan_context
    adata = ads.adata_dic[ads.active_id]
    func_kwargs = request.model_dump(exclude_none=True)
    pval = func_kwargs.pop("specificity_cutoff", 0.05)
    res_key = func_kwargs.get("uns_key", "liana_res")
    
    pval_col = adata.uns[res_key].columns[-1]
    func_kwargs["filter_fun"] = lambda x: x[pval_col] <= pval
    
    func_kwargs = {k: v for k, v in func_kwargs.items() 
                  if k in inspect.signature(li.pl.circle_plot).parameters}
    ax = li.pl.circle_plot(adata, **func_kwargs)
    func_kwargs["filter_fun"] = f"{pval}"
    fig_path = savefig(ax, "ccc_circle", **func_kwargs)
    add_op_log(adata, li.pl.circle_plot, func_kwargs)
    return {"figpath": str(fig_path)}


@pl_mcp.tool()
async def dotplot(request: DotPlotModel, ctx: Context):
    """Visualize cell-cell communication interactions using a dotplot."""
    result = await forward_request("ccc_dotplot", request.model_dump())
    if result is not None:
        return result
    
    ads = ctx.request_context.lifespan_context
    adata = ads.adata_dic[ads.active_id]

    func_kwargs = filter_args(request, li.pl.dotplot)
    pval = func_kwargs.pop("specificity_cutoff", 0.05)
    res_key = func_kwargs.get("uns_key", "liana_res")
    pval_col = adata.uns[res_key].columns[-1]
    func_kwargs["filter_fun"] = lambda x: x[pval_col] <= pval
    if func_kwargs.get("colour", None) is None:
        func_kwargs["colour"] = adata.uns[res_key].columns[-2]
    if func_kwargs.get("size", None) is None:
        func_kwargs["size"] = adata.uns[res_key].columns[-1]
    
    fig = li.pl.dotplot(adata, **func_kwargs)
    func_kwargs["filter_fun"] = f"{pval}"
    fig_path = savefig(fig, "ccc_dotplot", **func_kwargs)
    func_kwargs["filter_fun"] = f"lambda x: x[{pval_col}] <= {pval}"
    add_op_log(adata, li.pl.dotplot, func_kwargs)
    return {"figpath": str(fig_path)}
