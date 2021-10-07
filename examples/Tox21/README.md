# Tox21 Descriptors

Download the Tox21 datasets (train, val, test) from the
 [Tox21 Data Challenge Website](https://tripod.nih.gov/tox21/challenge/data.jsp)

## Datasets

```yaml
train: tox21_10k_data_all.sdf
val:   tox21_10k_challenge_score.sdf
test:  tox21_10k_challenge_test.smiles
```

## Usage

```bash
compute_descriptors \
    --config-dir $(pwd) \
    --config-name tox21 \
    data.dataset_path=$(pwd)/tox21_10k_data_all.sdf \
    ray.dashboard_port=8265
```

## Remarks

Always check whether all compounds could be parsed by RDKit.
It's possible to reparse (read and export again) RDKit-unparsable SMILES using OpenBabel.
The new OpenBabel Exported SMILES are then often readable by RDKit.
