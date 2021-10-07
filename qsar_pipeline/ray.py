import os
from typing import List

import pandas as pd
import ray
from omegaconf import OmegaConf

from qsar_pipeline import DescriptorCalculatorActor, utils


def compute_descriptors(molitem_list: List[utils.MolItem], cfg: OmegaConf):
    assert not OmegaConf.is_missing(cfg, "ray.dashboard_port")

    try:
        username = os.environ["USER"]
        ray.init(
            _temp_dir=f"/tmp/{username}/ray",
            dashboard_port=cfg.ray.dashboard_port,
            log_to_driver=False,
        )

        actors = [
            DescriptorCalculatorActor.remote(cfg.data)  # pylint: disable=no-member
            for _ in range(cfg.ray.num_actors)
        ]
        actor_pool = ray.util.ActorPool(actors)

        result_list = list(
            actor_pool.map_unordered(
                lambda actor, molitem: actor.calculate_descriptors.remote(
                    molitem,
                    standardize_mol=cfg.data.standardize_mol,
                    return_smiles=cfg.data.return_smiles,
                    return_properties=cfg.data.return_properties,
                ),
                molitem_list,
            )
        )
        df_results = pd.DataFrame.from_records(result_list, index="index")
        return df_results
    except Exception as exception:
        raise RuntimeError("Exception encountered") from exception
    finally:
        ray.shutdown()
