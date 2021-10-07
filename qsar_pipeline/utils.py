import gzip
import pickle
import string
from pathlib import Path
from typing import List, NamedTuple

import numpy as np
import pandas as pd
from omegaconf.omegaconf import OmegaConf
from openbabel import openbabel
from rdkit import Chem

from qsar_pipeline.calculators.rdkit_calculator import (
    morgan_fc_features,
    morgan_fp_features,
)


class MolItem(NamedTuple):
    index: int
    mol: Chem.Mol

    @classmethod
    def sanitize(cls, index: int, mol: Chem.Mol):
        mol = check_pickle5(mol)
        return cls(index, mol)


def check_pickle5(mol: Chem.Mol):
    if mol is not None:
        try:
            mol_pkl = pickle.dumps(mol, protocol=5)
            _ = pickle.loads(mol_pkl)
            return mol
        except RuntimeError:
            return None
    return None


def smiles_to_mol(smiles: str):
    if smiles is not None and smiles is not np.nan:
        if isinstance(smiles, (list, np.ndarray)):
            smiles = smiles[0]  # if multiple SMILES available use first one
        mol = Chem.MolFromSmiles(smiles)
        return mol
    return None


def _convert_smiles(smiles, obconversion):
    if isinstance(smiles, str):
        mol = openbabel.OBMol()
        obconversion.ReadString(mol, smiles)
        smiles_reparsed = obconversion.WriteString(mol)
        smiles_reparsed = smiles_reparsed.strip()

        # try to re-read smiles using RDKit
        _ = Chem.MolFromSmiles(smiles_reparsed)
        return smiles_reparsed
    return None


def reparse_smiles(smiles_list):
    obconversion = openbabel.OBConversion()
    obconversion.SetInAndOutFormats("smi", "smi")

    smiles_list_reparsed = []

    for idx, smiles_sublist in smiles_list:
        if isinstance(smiles_sublist, list):
            smiles_sublist_reparsed = []
            for smiles in smiles_sublist:
                smiles_reparsed = _convert_smiles(smiles, obconversion)
                smiles_sublist_reparsed.append(smiles_reparsed)

            smiles_list_reparsed.append((idx, smiles_sublist_reparsed))
        else:
            # smiles_sublist is single SMILE
            smiles = smiles_sublist
            smiles_reparsed = _convert_smiles(smiles, obconversion)
            smiles_list_reparsed.append((idx, smiles_reparsed))

    return smiles_list_reparsed


def parse_smiles(
    smi_path: Path,
    smiles_column: str = "smiles",
    filetype: str = "csv",
    reparse: bool = False,
):
    smi_path = smi_path.resolve() if smi_path.is_symlink() else smi_path

    if filetype == "parquet":
        df_raw = pd.read_parquet(smi_path)
    elif filetype == "csv":
        df_raw = pd.read_csv(smi_path)
    elif filetype == "tsv":
        df_raw = pd.read_csv(smi_path, sep="\t")

    smiles_list = list(df_raw[smiles_column].items())

    if reparse is True:
        smiles_list = reparse_smiles(smiles_list)

    mol_list = [
        MolItem.sanitize(index, smiles_to_mol(smiles)) for index, smiles in smiles_list
    ]

    return mol_list, df_raw


def parse_sdf(sdf_path: Path, compression: str = "none"):
    sdf_path = sdf_path.resolve() if sdf_path.is_symlink() else sdf_path

    if compression == "none":
        suppl = Chem.SDMolSupplier(sdf_path.as_posix())
        for index, mol in enumerate(suppl):
            yield MolItem.sanitize(index, mol)
    elif compression == "gzip":
        with gzip.open(sdf_path) as sdf:
            suppl = Chem.ForwardSDMolSupplier(sdf)
            for index, mol in enumerate(suppl):
                yield MolItem.sanitize(index, mol)


def parse_sdf_as_dataframe(sdf_path: Path, compression: str = "none"):
    df_raw = [
        {"index": index, "molblock": Chem.MolToMolBlock(mol), **mol.GetPropsAsDict()}
        for index, mol in parse_sdf(sdf_path, compression)
    ]
    df_raw = pd.DataFrame.from_records(df_raw, index="index")
    return df_raw


def _check_argument(feature_config: OmegaConf, config_option: str, feature_name: str):
    option_value = OmegaConf.select(
        feature_config, config_option, throw_on_missing=True
    )
    if option_value is None:
        raise ValueError(
            (
                f"Mandatory option `{config_option}` "
                f"in feature definition `{feature_name}` is not set"
            )
        )


def check_missing(cfg_features: OmegaConf):
    for feature_name, feature_config in cfg_features.items():
        if feature_config.type in morgan_fp_features:
            _check_argument(feature_config, "radius", feature_name)
            _check_argument(feature_config, "nBits", feature_name)

        if feature_config.type in morgan_fc_features:
            _check_argument(feature_config, "radius", feature_name)


def get_filetype_suffix(filetype_suffixes: List[str]):
    filetype_suffix_set = {".smi", ".csv", ".tsv"} & set(filetype_suffixes)
    if len(filetype_suffix_set) > 1:
        raise ValueError("Conflicting filetypes. Cannot determine filetype of input")

    filetype_suffix = filetype_suffix_set.pop()
    filetype_suffix = filetype_suffix.translate(
        str.maketrans("", "", string.punctuation)
    )

    return filetype_suffix
