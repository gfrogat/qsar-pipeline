import mordred
from mordred import descriptors as mordred_descriptors
from omegaconf.omegaconf import OmegaConf


def _compute_mordred_descriptors(feature_name, calculator, mol):
    mordred_results = calculator(mol).fill_missing()
    mordred_results = {
        f"{feature_name}::mordred_descriptors::{key}": value
        for key, value in mordred_results.items()
    }
    return mordred_results


class MordredCalculator:
    _supported_features = ["mordred", "mordred_descriptors"]

    @staticmethod
    def supported_features(return_set=False):
        features = MordredCalculator._supported_features
        features = set(features) if return_set else features
        return features

    def __init__(self, cfg_features: OmegaConf):
        self.calculators = {}

        for feature_name, feature_config in cfg_features.items():
            if feature_config.type in MordredCalculator.supported_features():
                self.calculators[feature_name] = mordred.Calculator(
                    mordred_descriptors, ignore_3D=True, version="1.2.0"
                )

    def __call__(self, mol):
        mordred_results = {}

        for feature_name, calculator in self.calculators.items():
            mordred_results.update(
                _compute_mordred_descriptors(feature_name, calculator, mol)
            )

        return mordred_results
