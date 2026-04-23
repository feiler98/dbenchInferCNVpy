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
    return {p.stem: p for p in target_path.rglob("*__hg_38.txt")}


def csvs_to_adatas(target_path: Path) -> dict:
    """
    Generates a dictionary with adata and their respective reference catalogue of normal cells (cell_names).
    """
    dict_hg38_desc = get_hg_38_desc_paths(target_path)
    dict_accepted_files = {}
    for k, path_txt in dict_hg38_desc.items():
        path_rcm = Path(target_path) / f"{k}__RCM.csv"
        if path_rcm.exists():
            adata = sc.read_csv(path_rcm).T
            adata.obs["cell_names"] = adata.obs.index
            with open(path_txt, "r") as f:
                list_norm_cells = list(map(lambda x: x.replace("\n", ""), f.readlines()))
                dict_accepted_files[k] = {"adata": adata,
                                          "reference_key": "cell_names",
                                          "reference_cat":list_norm_cells}
    return dict_accepted_files


def run_py_infercnv(path_target: Path, path_out_data: Path, kwargs: dict = {}) -> None:
    """
    Main function for running infercnvpy for benchmarking.

    Parameters
    ----------
    path_target: Path
        Directory with all datasets for benchmarking.
    path_out_data: Path
        Directory where to save the results for benchmarking.
    kwargs: dict
        Key-word-arguments (= kwargs) for infercnvpy.tl.infercnv() function.

    Returns
    -------
    pd.DataFrame
        Returns inferred copy number variations as table.
    """
    dict_files = csvs_to_adatas(path_target)
    for tag_dataset, dict_data in dict_files.items():
        str_kwargs = ";".join([f"{list(x)[0]},{y}" for x, y in kwargs.items()])
        file_name = f"{tag_dataset}__{str_kwargs}__infercnvpy"
        data_save_path = path_out_data / file_name
        data_save_path.mkdir(exist_ok=True)

        adata = dict_data["adata"]
        # gencode hg38 file needed for providing "start", "end", & "chr" information
        cnv.io.genomic_position_from_gtf("gencode.v38.annotation.gtf", adata)

        @benchmark_method(data_save_path)
        def run_infercnvpy_func(adata, dict_data, kwargs):
            cnv.tl.infercnv(adata, calculate_gene_values=True,
                            reference_key=dict_data["reference_key"],
                            reference_cat=dict_data["reference_cat"],
                            **kwargs)

        run_infercnvpy_func(adata, dict_data, kwargs)
        cnv_idx = list(adata.obs.index)

        df_csv_pre = pd.DataFrame(data=adata.layers["gene_values_cnv"], index=cnv_idx).T
        df_csv = pd.concat([adata.var.reset_index(), df_csv_pre], axis=1).dropna().drop("index", axis=1).set_index("gene_name")
        df_csv.to_csv(data_save_path / f"{file_name}.csv")


if __name__ == "__main__":

    # matrix of possible infercnv_py hyperparameter kwargs
    kwargs_gridsearch = {
        "n_jobs": [5, 10, 20],
        "step": [1, 5, 10, 20],
        "window_size": [10, 25, 50, 100, 200, 500],
        "dynamic_threshold": [None, 1, 1.5, 2, 3],
        "chunksize": [10, 50, 100, 500, 1000],
    }

    path_in, path_out = val_build_project()
    run_py_infercnv(path_in, path_out, kwargs={"n_jobs": 1, "chunksize": 100})  # 1 core computing standard params adjusted for calculating gene values
    list_kwargs = grid_by_dict(kwargs_gridsearch)
    for kwarg_opt in list_kwargs:
        print(f"InferCNVpy running with hyperparameters: {kwarg_opt}")
        run_py_infercnv(path_in, path_out, kwargs=kwarg_opt)