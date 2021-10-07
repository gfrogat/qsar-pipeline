import logging
from pathlib import Path

import hydra
from omegaconf import OmegaConf

import qsar_pipeline.ray as qsar_pipeline_ray
import qsar_pipeline.utils as utils


@hydra.main(config_path="../config", config_name="config.yaml")
def compute_descriptors(cfg: OmegaConf):
    logger = logging.getLogger("QSAR-Pipeline")

    dataset_path = Path(cfg.data.dataset_path)
    dataset_path = dataset_path.resolve() if dataset_path.is_symlink() else dataset_path
    filetype_suffixes = dataset_path.suffixes

    if ".sdf" in filetype_suffixes:
        logger.info("Loading molecules from SDF")
        molitem_list = list(
            utils.parse_sdf(dataset_path, cfg.data.sdf_dataset_compression)
        )
        df_raw = None
    elif {".smi", ".csv", ".tsv"} & set(filetype_suffixes):
        filetype_suffix = utils.get_filetype_suffix(filetype_suffixes)
        logger.info("Loading molecules from %s", filetype_suffix)
        molitem_list, df_raw = utils.parse_smiles(
            dataset_path,
            cfg.data.smiles_column,
            filetype=filetype_suffix,
            reparse=cfg.data.reparse_smiles,
        )
    elif ".parquet" in filetype_suffixes:
        logger.info("Loading molecules from Parquet")
        molitem_list, df_raw = utils.parse_smiles(
            dataset_path, cfg.data.smiles_column, filetype="parquet"
        )
    else:
        raise ValueError("Unknown filetype")

    logger.info("Checking config")
    utils.check_missing(cfg.data.chem_features)

    logger.info("Computing descriptors")
    df_results = qsar_pipeline_ray.compute_descriptors(molitem_list, cfg)

    if df_raw is not None:
        df_raw.columns = list(
            map(lambda colname: "".join(["raw::", colname]), df_raw.columns)
        )
        df_results = df_results.join(df_raw, on="index")

    logger.info("Storing results")
    df_results.to_parquet(Path(cfg.data.result_path).resolve())

    logger.info("Storing config")
    OmegaConf.save(cfg, f=Path("pipeline_config.yaml"))


def entry():
    compute_descriptors()  # pylint: disable=no-value-for-parameter


if __name__ == "__main__":
    compute_descriptors()  # pylint: disable=no-value-for-parameter
