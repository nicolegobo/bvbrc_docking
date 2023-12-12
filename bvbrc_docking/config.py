from typing import Literal, Optional, Union

from pydantic import Field

from bvbrc_docking.utils import BaseModel


class fredConfig(BaseModel):
    name: Literal["fred"] = "fred"
    receptor_pdb: Optional[str] = Field(
        ..., description="pdb path for the input protein"
    )
    drug_dbs: str = Field(..., description="smi path for the ligands")
    n_cpus: int = Field(..., description="number of MPI threads")
    fred_path: Optional[str]


class fredPartialConfig(BaseModel):
    """used for colmena runs"""

    name: Literal["fred_partial"] = "fred_partial"
    drug_dbs: str = Field(..., description="smi path for the ligands")
    n_cpus: int = Field(..., description="number of MPI threads")
    fred_path: Optional[str] = Field(..., description="path to fred exec")


class DiffDockConfig(BaseModel):
    name: Literal["diffdock"] = "diffdock"
    receptor_pdb: Optional[str] = Field(
        ..., description="pdb path for the input protein"
    )
    drug_dbs: str = Field(..., description="smi path for the ligands")
    work_dir: str = Field(..., description="installed path for diffdock")
    output_path: str = Field(..., description="output path for the docking result")
    top_n: Optional[int] = Field(
        ...,
        description="number of top N candidates for each protein-ligand pair",
    )


class DiffDockPartialConfig(BaseModel):
    """used for colmena runs"""

    name: Literal["diffdock_partial"] = "diffdock_partial"
    drug_dbs: str = Field(..., description="smi path for the ligands")
    work_dir: str = Field(..., description="installed path for diffdock")
    output_path: str = Field(..., description="output path for the docking result")
    top_n: Optional[int] = Field(
        ...,
        description="number of top N candidates for each protein-ligand pair",
    )


DockConfig = Union[fredConfig, DiffDockConfig, fredPartialConfig, DiffDockPartialConfig]
