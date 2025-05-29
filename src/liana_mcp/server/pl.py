from fastmcp import FastMCP
from fastmcp.exceptions import ToolError
import liana as li
from ..schema.pl import *
from scmcp_shared.util import add_op_log, filter_args, forward_request, savefig, get_ads
from scmcp_shared.logging_config import setup_logger
from scmcp_shared.schema import AdataInfo
logger = setup_logger()


pl_mcp = FastMCP("lianaMCP-pl-Server")


@pl_mcp.tool()
def circle_plot(
    request: CirclePlotModel, 
    adinfo: AdataInfo = AdataInfo()
):
    """Visualize cell-cell communication network using a circular plot."""
    try:
        result = forward_request("pl_circle_plot", request, adinfo)
        if result is not None:
            return result
        ads = get_ads()
        adata = ads.get_adata(adinfo=adinfo)
        func_kwargs = filter_args(request, li.pl.circle_plot)
        pval = func_kwargs.pop("specificity_cutoff", 0.05)
        res_key = func_kwargs.get("uns_key", "liana_res")
        pval_col = adata.uns[res_key].columns[-1]
        if "filter_fun" in func_kwargs:
            func_kwargs["filter_fun"] = lambda x: x[pval_col] <= pval
        ax = li.pl.circle_plot(adata, **func_kwargs)
        func_kwargs["filter_fun"] = f"{pval}"
        fig_path = savefig(ax, li.pl.circle_plot, **func_kwargs)
        add_op_log(adata, li.pl.circle_plot, func_kwargs, adinfo)
        return {"figpath": str(fig_path)}
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)


@pl_mcp.tool()
def ccc_dotplot(
    request: DotPlotModel, 
    adinfo: AdataInfo = AdataInfo()
):
    """Visualize cell-cell communication interactions using a dotplot."""
    try:
        result = forward_request("ccc_dotplot", request, adinfo)
        if result is not None:
            return result
        
        ads = get_ads()
        adata = ads.get_adata(adinfo=adinfo)
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
        fig_path = savefig(fig, li.pl.dotplot, **func_kwargs)
        func_kwargs["filter_fun"] = f"lambda x: x[{pval_col}] <= {pval}"
        add_op_log(adata, li.pl.dotplot, func_kwargs, adinfo)
        return {"figpath": str(fig_path)}
    except ToolError as e:
        raise ToolError(e)
    except Exception as e:
        if hasattr(e, '__context__') and e.__context__:
            raise ToolError(e.__context__)
        else:
            raise ToolError(e)
