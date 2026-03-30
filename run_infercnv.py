# imports
# ----------------------------------------------------------------------------------------------------------------------
import pandas as pd
from pathlib import Path
from pyomics.utils import benchmark_method
import infercnvpy as cnv
import scanpy as sc
import itertools
# ----------------------------------------------------------------------------------------------------------------------

def grid_by_dict(pars_dict: dict) -> list:
    keys=pars_dict.keys()
    combinations = itertools.product(*pars_dict.values())
    list_of_kwargs = [dict(zip(keys, cc)) for cc in combinations]
    return list_of_kwargs


def val_build_project() -> (Path, Path):
    cwd_path = Path.cwd()
    print(f"Current working directory of running script {Path(__file__).name}: {cwd_path}")
    path_out = cwd_path / "app" / "out"
    path_in = cwd_path / "data_input"

    if not path_in.exists():
        raise ValueError(f"Data dir '{str(path_in)}' does not exist!")

    if not path_out.exists():
        path_out.mkdir(parents=True, exist_ok=True)
        print(f"Data out-dir '{str(path_out)}' has been created...")

    return path_in, path_out


def get_hg_38_desc_paths(target_path: Path) -> dict:
    """
    These fetched .txt files correlate to .csv RCM files --> describe normal cells within the datasets.
    """
    return {p.stem: p for p in target_path.rglob("*__hg_38__RCM.txt")}


def csvs_to_adatas(target_path: Path) -> dict:
    """
    Generates a dictionary with adata and their respective reference catalogue of normal cells (cell_names).
    """
    dict_hg38_desc = get_hg_38_desc_paths(target_path)
    dict_accepted_files = {}
    for k, path_txt in dict_hg38_desc.items():
        path_rcm = Path(target_path) / f"{k}.csv"
        if path_rcm.exists():
            adata = sc.read_csv(path_rcm).T
            adata.obs["cell_names"] = adata.obs.index
            with open(path_txt, "r") as f:
                list_norm_cells = list(map(lambda x: x.replace("\n", ""), f.readlines()))
                dict_accepted_files[k] = {"adata": adata,
                                          "reference_key": "cell_names",
                                          "reference_cat":list_norm_cells}
    return dict_accepted_files


def run_py_infercnv(path_target: Path, kwargs: dict = dict) -> pd.DataFrame:
    """
    Main function for running infercnvpy for benchmarking.

    Parameters
    ----------
    path_target: Path
        Directory with all datasets for benchmarking.
    kwargs: dict
        Key-word-arguments (= kwargs) for infercnvpy.tl.infercnv() function.

    Returns
    -------
    pd.DataFrame
        Returns inferred copy number variations as table.
    """
    adata = ""

    cnv.tl.infercnv(adata, **kwargs)
    cnv_X = adata.obsm["X_cnv"].toarray()

    # shape of the numpy-array cnv_X
    shape_tuple = cnv_X.shape
    cnv_idx = list(adata.obs.index)

    # col_tags
    chr_start_pos_list = list(adata.uns["cnv"]["chr_pos"].values())
    chr_tags = list(adata.uns["cnv"]["chr_pos"].keys())
    chr_start_pos_list_pop_first = chr_start_pos_list[1::]
    chr_start_pos_list_pop_first.append(shape_tuple[1])
    dict_chr_end_start_pos = {chr_tag :range(int(e ) -int(s)) for chr_tag, s, e in zip(chr_tags, chr_start_pos_list, chr_start_pos_list_pop_first)}

    cnv_col = []
    for chr_tag, range_list in dict_chr_end_start_pos.items():
        for pos in range_list:
            cnv_col.append(f"{chr_tag}__rel_pos_{pos}")
    return pd.DataFrame(data=cnv_X, columns=cnv_col, index=cnv_idx).T


if __name__ == "__main__":
    pass