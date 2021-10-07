import ray
from omegaconf import OmegaConf
from rdkit import Chem
from rdkit.Chem import MolStandardize

from qsar_pipeline import utils
from qsar_pipeline.calculators import (
    MordredCalculator,
    RDKitCalculator,
    ToxicophoresCalculator,
)

# required for ray
Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)


@ray.remote(num_cpus=1)
class DescriptorCalculatorActor:
    def __init__(self, cfg: OmegaConf):
        self.cfg = cfg
        self.descriptor_calculator = DescriptorCalculator(cfg)

    def get_cfg(self):
        return self.cfg

    @ray.method(num_returns=1)
    def calculate_descriptors(
        self,
        molitem: utils.MolItem,
        standardize_mol: bool = True,
        return_smiles: bool = False,
        return_properties: bool = False,
    ):
        return self.descriptor_calculator(
            molitem, standardize_mol, return_smiles, return_properties
        )


class DescriptorCalculator:
    def __init__(self, cfg: OmegaConf):
        self.rdkit_calculator = DescriptorCalculator.setup_rdkit_calculator(cfg)
        self.mordred_calculator = DescriptorCalculator.setup_mordred_calculator(cfg)
        self.toxicophores_calculator = (
            DescriptorCalculator.setup_toxicophores_calculator(cfg)
        )

    @staticmethod
    def setup_rdkit_calculator(cfg: OmegaConf):
        chem_features = {value.type for value in cfg.chem_features.values()}
        rdkit_features = RDKitCalculator.supported_features(return_set=True)

        if len(chem_features & rdkit_features) > 0:
            return RDKitCalculator(cfg.chem_features)

        return None

    @staticmethod
    def setup_mordred_calculator(cfg: OmegaConf):
        chem_features = {value.type for value in cfg.chem_features.values()}
        mordred_features = MordredCalculator.supported_features(return_set=True)

        if len(chem_features & mordred_features) > 0:
            return MordredCalculator(cfg.chem_features)

        return None

    @staticmethod
    def setup_toxicophores_calculator(cfg: OmegaConf):
        chem_features = {value.type for value in cfg.chem_features.values()}
        toxicophores_features = ToxicophoresCalculator.supported_features(
            return_set=True
        )

        if len(chem_features & toxicophores_features) > 0:
            return ToxicophoresCalculator(cfg.chem_features)

        return None

    def __call__(
        self,
        molitem: utils.MolItem,
        standardize_mol: bool = True,
        return_smiles: bool = True,
        return_properties: bool = False,
    ):
        results = {"index": molitem.index}
        mol = molitem.mol

        if mol is None:
            return results

        if standardize_mol:
            mol = MolStandardize.rdMolStandardize.Cleanup(mol)

        if return_smiles:
            standardized_str = "_std" if standardize_mol else ""
            results[f"smiles{standardized_str}"] = Chem.MolToSmiles(mol)

        if self.rdkit_calculator is not None:
            rdkit_results = self.rdkit_calculator(mol)
            results.update(rdkit_results)

        if self.mordred_calculator is not None:
            mordred_results = self.mordred_calculator(mol)
            results.update(mordred_results)

        if self.toxicophores_calculator is not None:
            toxicophores_results = self.toxicophores_calculator(mol)
            results.update(toxicophores_results)

        if return_properties:
            properties = {
                f"properties::{prop_name}": prop
                for prop_name, prop in mol.GetPropsAsDict().items()
            }
            results.update(properties)

        return results
