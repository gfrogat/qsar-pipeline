from functools import partial

import numpy as np
from omegaconf import OmegaConf
from rdkit.Chem import rdMolDescriptors


def _compute_maccs_keys(mol, feature_name, feature_config):
    kwargs = OmegaConf.select(feature_config, "kwargs", default={})
    maccs_keys = rdMolDescriptors.GetMACCSKeysFingerprint(mol, **kwargs)
    maccs_keys = np.array(maccs_keys.GetOnBits())
    return {f"{feature_name}::maccs_keys@167": maccs_keys}


def _compute_morgan_fp(mol, feature_name, feature_config):
    radius = OmegaConf.select(feature_config, "radius", throw_on_missing=True)
    nbits = OmegaConf.select(feature_config, "nBits", throw_on_missing=True)
    kwargs = OmegaConf.select(feature_config, "kwargs", default={})

    morgan_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol, radius=radius, nBits=nbits, **kwargs
    )
    morgan_fp = np.array(morgan_fp.GetOnBits())

    return {f"{feature_name}::morgan_fp{radius*2}@{nbits}": morgan_fp}


def _compute_morgan_fc(mol, feature_name, feature_config):
    radius = OmegaConf.select(feature_config, "radius", throw_on_missing=True)
    kwargs = OmegaConf.select(feature_config, "kwargs", default={})

    morgan_fc = rdMolDescriptors.GetMorganFingerprint(
        mol, radius=radius, **kwargs
    ).GetNonzeroElements()

    return {
        f"{feature_name}::morgan_fc{radius*2}_keys": list(morgan_fc.keys()),
        f"{feature_name}::morgan_fc{radius*2}_values": list(morgan_fc.values()),
    }


maccs_features = {"maccs": _compute_maccs_keys, "maccs_keys": _compute_maccs_keys}

morgan_fp_features = {
    "morgan_fp": _compute_morgan_fp,
    "ecfp": _compute_morgan_fp,
}
morgan_fc_features = {
    "morgan_fc": _compute_morgan_fc,
    "ecfc": _compute_morgan_fc,
}


class RDKitCalculator:
    _supported_features = {**maccs_features, **morgan_fp_features, **morgan_fc_features}

    @staticmethod
    def supported_features(return_set=False):
        features = list(RDKitCalculator._supported_features.keys())
        features = set(features) if return_set else features
        return features

    def __init__(self, cfg_features: OmegaConf):
        self.transforms = []
        self.cfg_features = cfg_features
        for feature_name, feature_config in cfg_features.items():
            if feature_config.type in RDKitCalculator.supported_features():
                self.transforms.append(
                    partial(
                        RDKitCalculator._supported_features[feature_config.type],
                        feature_name=feature_name,
                        feature_config=feature_config,
                    )
                )

    def __call__(self, mol):
        results = {}

        for transform in self.transforms:
            results.update(transform(mol))

        return results

    def __repr__(self):
        return OmegaConf.to_yaml(self.cfg_features)
