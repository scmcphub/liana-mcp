from fastmcp import FastMCP, Context
from fastmcp.exceptions import ToolError
import liana as li
import inspect
from pathlib import Path
import os
from ..schema.pl import *
from scmcp_shared.util import add_op_log, filter_args, forward_request, set_fig_path, get_ads
from scmcp_shared.logging_config import setup_logger

logger = setup_logger()


pl_mcp = FastMCP("lianaMCP-pl-Server")


@pl_mcp.tool()
async def circle_plot(
    request: CirclePlotModel, 
):
    """Visualize cell-cell communication network using a circular plot."""
    try:
        result = await forward_request("pl_circle_plot", request)
        if result is not None:
            return result
        ads = get_ads()
        adata = ads.get_adata(request=request)
        func_kwargs = filter_args(request, li.pl.circle_plot)
        pval = func_kwargs.pop("specificity_cutoff", 0.05)
        res_key = func_kwargs.get("uns_key", "liana_res")
        pval_col = adata.uns[res_key].columns[-1]
        if "filter_fun" in func_kwargs:
            func_kwargs["filter_fun"] = lambda x: x[pval_col] <= pval
        ax = li.pl.circle_plot(adata, **func_kwargs)
        func_kwargs["filter_fun"] = f"{pval}"
        fig_path = set_fig_path(ax, li.pl.circle_plot, **func_kwargs)
        add_op_log(adata, li.pl.circle_plot, func_kwargs)
        return {"figpath": str(fig_path)}
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)


@pl_mcp.tool()
async def ccc_dotplot(
    request: DotPlotModel, 
):
    """Visualize cell-cell communication interactions using a dotplot."""
    try:
        result = await forward_request("ccc_dotplot", request)
        if result is not None:
            return result
        
        ads = get_ads()
        adata = ads.get_adata(request=request)
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
        fig_path = set_fig_path(fig, li.pl.dotplot, **func_kwargs)
        func_kwargs["filter_fun"] = f"lambda x: x[{pval_col}] <= {pval}"
        add_op_log(adata, li.pl.dotplot, func_kwargs)
        return {"figpath": str(fig_path)}
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)
