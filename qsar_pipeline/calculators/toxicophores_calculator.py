import numpy as np
from omegaconf.omegaconf import OmegaConf
from rdkit import Chem

from qsar_pipeline.tox_smarts import tox_smarts


def _init_patterns():
    smarts_patts = []
    for smarts in tox_smarts:
        smarts_query = smarts.split(" ")
        query_components = []

        for element in smarts_query:
            if element in ["NOT", "AND", "OR"]:
                query_components.append((element.lower(), False))
            else:
                query_components.append((Chem.MolFromSmarts(element), True))

        smarts_patts.append(query_components)

    return smarts_patts


def _compute_tox_keys(feature_name, smarts_patterns, mol):
    tox_keys = []
    nbits = len(smarts_patterns)

    for idx, query_components in enumerate(smarts_patterns):
        query = []

        for (patt_op, is_patt) in query_components:
            if is_patt is True:
                patt = patt_op
                is_match = mol.HasSubstructMatch(patt)
                query.append(str(is_match))
            else:
                op = patt_op  # pylint: disable=invalid-name
                query.append(op)

        query = " ".join(query)
        # eval is only used for evaluating the boolean expression
        query_result = eval(query)  # pylint: disable=eval-used

        if query_result is True:
            tox_keys.append(idx)

    return {f"{feature_name}::tox_keys@{nbits}": np.array(tox_keys)}


class ToxicophoresCalculator:
    _supported_features = ["tox", "tox_keys", "toxicophores"]

    @staticmethod
    def supported_features(return_set=False):
        features = ToxicophoresCalculator._supported_features
        features = set(features) if return_set else features
        return features

    def __init__(self, cfg_features: OmegaConf):
        self.smarts_patterns = {}
        for feature_name, feature_config in cfg_features.items():
            if feature_config.type in ToxicophoresCalculator.supported_features():
                self.smarts_patterns[feature_name] = _init_patterns()

    def __call__(self, mol):
        toxicophores_results = {}

        for feature_name, smarts_patterns in self.smarts_patterns.items():
            toxicophores_results.update(
                _compute_tox_keys(feature_name, smarts_patterns, mol)
            )
        return toxicophores_results
