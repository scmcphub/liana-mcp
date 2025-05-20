from pydantic import (
    Field,
    ValidationInfo,
    computed_field,
    field_validator,
    model_validator,
)
from typing import Optional, List, Dict, Union, Literal, Tuple
from scmcp_shared.schema import AdataModel

class ListCCCMethodModel(AdataModel):
    """ListCCCMethodModel"""    
    pass    


class RankAggregateModel(AdataModel):
    """Input schema for LIANA's rank_aggregate method for cell-cell communication analysis."""
    
    groupby: str = Field(
        ...,  # Required field
        description="Key to be used for grouping or clustering cells (e.g., cell type annotations)."
    )
    
    resource_name: str = Field(
        default="consensus",
        description="Name of the resource to be used for ligand-receptor inference. See `li.rs.show_resources()` for available resources."
    )
    
    expr_prop: float = Field(
        default=0.1,
        description="Minimum expression proportion for the ligands and receptors in the corresponding cell identities. Set to 0 to return unfiltered results."
    )
    
    min_cells: int = Field(
        default=5,
        description="Minimum cells per cell identity to be considered for downstream analysis."
    )
    
    base: float = Field(
        default=2.718281828459045,  # e
        description="Exponent base used to reverse the log-transformation of the matrix. Relevant only for the `logfc` method."
    )
    
    aggregate_method: Literal["rra", "mean"] = Field(
        default="rra",
        description="Method aggregation approach: 'mean' for mean rank, 'rra' for RobustRankAggregate."
    )
    
    return_all_lrs: bool = Field(
        default=False,
        description="Whether to return all ligand-receptor pairs, or only those that surpass the expr_prop threshold."
    )
    
    key_added: str = Field(
        default="liana_res",
        description="Key under which the results will be stored in adata.uns."
    )
    
    use_raw: Optional[bool] = Field(
        default=True,
        description="Use raw attribute of adata if present."
    )
    
    layer: Optional[str] = Field(
        default=None,
        description="Layer in AnnData.layers to use. If None, use AnnData.X."
    )
    
    de_method: str = Field(
        default="t-test",
        description="Differential expression method used to rank genes according to 1vsRest."
    )
    
    n_perms: int = Field(
        default=1000,
        description="Number of permutations for permutation-based methods. If None, no permutation testing is performed."
    )

    n_jobs: int = Field(
        default=1,
        description="Number of jobs to run in parallel."
    )


class CCCModel(AdataModel):
    """Input schema for LIANA's cell-cell communication analysis."""
    
    method: Literal[
        "singlecellsignalr", 
        "connectome", 
        "cellphonedb", 
        "natmi", 
        "logfc", 
        "cellchat", 
        "geometric_mean", 
        "scseqcomm"
    ] = Field(
        default="cellphonedb",
        description="cell-cell communication method"
    )
    
    groupby: str = Field(
        ...,  # Required field
        description="Key to be used for grouping cells (e.g., cell type annotations)."
    )
    
    resource_name: str = Field(
        default="consensus",
        description="Name of the resource to be used for ligand-receptor inference. See `li.rs.show_resources()` for available resources."
    )
    
    expr_prop: float = Field(
        default=0.1,
        description="Minimum expression proportion for the ligands and receptors in the corresponding cell identities. Set to 0 to return unfiltered results."
    )
    
    min_cells: int = Field(
        default=5,
        description="Minimum cells per cell identity to be considered for downstream analysis."
    )
    
    base: float = Field(
        default=2.718281828459045,  # e
        description="Exponent base used to reverse the log-transformation of the matrix. Relevant only for the `logfc` method."
    )
    
    return_all_lrs: bool = Field(
        default=False,
        description="Whether to return all ligand-receptor pairs, or only those that surpass the expr_prop threshold."
    )
    
    key_added: str = Field(
        default="liana_res",
        description="Key under which the results will be stored in adata.uns."
    )
    
    use_raw: Optional[bool] = Field(
        default=True,
        description="Use raw attribute of adata if present."
    )
    
    layer: Optional[str] = Field(
        default=None,
        description="Layer in AnnData.layers to use. If None, use AnnData.X."
    )
    
    de_method: str = Field(
        default="t-test",
        description="Differential expression method used to rank genes according to 1vsRest."
    )
    
    n_perms: int = Field(
        default=1000,
        description="Number of permutations for the permutation test. Relevant for CellPhoneDB method."
    )
    
    seed: int = Field(
        default=1337,
        description="Random seed for reproducibility."
    )
    n_jobs: int = Field(
        default=1,
        description="Number of jobs to run in parallel."
    )
    inplace: bool = Field(
        default=True,
        description="Whether to store results in place, or return them."
    )
    supp_columns: Optional[List[str]] = Field(
        default=None,
        description="Additional columns to be added from methods in liana, or columns from scanpy.tl.rank_genes_groups."
    )

