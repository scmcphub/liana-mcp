from fastmcp import FastMCP, Context
import liana as li
import inspect
from pathlib import Path
import os
from ..schema.pl import *
from scmcp_shared.util import add_op_log, filter_args, forward_request, set_fig_path
from scmcp_shared.logging_config import setup_logger

logger = setup_logger()


pl_mcp = FastMCP("lianaMCP-pl-Server")


@pl_mcp.tool()
async def circle_plot(
    request: CirclePlotModel, 
    ctx: Context,
    dtype: str = Field(default="exp", description="the datatype of anndata.X"),
    sampleid: str = Field(default=None, description="adata sampleid for plotting")
):
    """Visualize cell-cell communication network using a circular plot."""
    try:
        result = await forward_request("pl_circle_plot", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(sampleid=sampleid, dtype=dtype)
        func_kwargs = filter_args(request, li.pl.circle_plot)
        pval = func_kwargs.pop("specificity_cutoff", 0.05)
        res_key = func_kwargs.get("uns_key", "liana_res")
        pval_col = adata.uns[res_key].columns[-1]
        if "filter_fun" in func_kwargs:
            func_kwargs["filter_fun"] = lambda x: x[pval_col] <= pval
        ax = li.pl.circle_plot(adata, **func_kwargs)
        func_kwargs["filter_fun"] = f"{pval}"
        fig_path = set_fig_path("circle_plot", fig=ax, **func_kwargs)
        add_op_log(adata, li.pl.circle_plot, func_kwargs)
        return {"figpath": str(fig_path)}
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e

@pl_mcp.tool()
async def ccc_dotplot(
    request: DotPlotModel, 
    ctx: Context,
    dtype: str = Field(default="exp", description="the datatype of anndata.X"),
    sampleid: str = Field(default=None, description="adata sampleid for plotting")
):
    """Visualize cell-cell communication interactions using a dotplot."""
    try:
        result = await forward_request("ccc_dotplot", request, sampleid=sampleid, dtype=dtype)
        if result is not None:
            return result
        
        ads = ctx.request_context.lifespan_context
        adata = ads.get_adata(sampleid=sampleid, dtype=dtype)

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
        fig_path = set_fig_path("dotplot", fig=fig, **func_kwargs)
        func_kwargs["filter_fun"] = f"lambda x: x[{pval_col}] <= {pval}"
        add_op_log(adata, li.pl.dotplot, func_kwargs)
        return {"figpath": str(fig_path)}
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise Exception(f"{str(e.__context__)}")
        else:
            raise e