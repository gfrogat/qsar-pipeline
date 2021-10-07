# Tox21Extended Descriptors

Download the file `tox21_10k_library_info.tsv.zip` from the
 [Tox21 Public Data](https://tripod.nih.gov/tox21/assays/) website.

## Usage

```bash
compute_descriptors \
    --config-dir $(pwd) \
    --config-name tox21_extended \
    data.dataset_path=$(pwd)/tox21_10k_library_info.tsv.zip \
    ray.dashboard_port=8265
```

## Remarks

Always check whether all compounds could be parsed by RDKit.
It's possible to reparse (read and export again) RDKit-unparsable SMILES using OpenBabel.
The new OpenBabel Exported SMILES are then often readable by RDKit.
