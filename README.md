# Small molecule standardization

The model outputs molecular representations at increasing levels of chemical abstraction, ranging from the full canonical structure to a ChEMBL-standardized parent, a stereochemistry-free form, the core Murcko scaffold, and a generic scaffold that removes atom-specific identities.

This model was incorporated on 2026-01-08.Last packaged on 2026-01-18.

## Information
### Identifiers
- **Ersilia Identifier:** `eos4k4f`
- **Slug:** `standardization`

### Domain
- **Task:** `Annotation`
- **Subtask:** `Property calculation or prediction`
- **Biomedical Area:** `Any`
- **Target Organism:** `Any`
- **Tags:** `Compound generation`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `5`
- **Output Consistency:** `Fixed`
- **Interpretation:** The model outputs representations at multiple chemical abstraction levels for a given molecule.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| canonical_smiles | string |  | Canonical SMILES for the molecule |
| standardized_smiles | string |  | ChEMBL-standardized parent SMILES |
| flattened_smiles | string |  | Standardized SMILES without stereochemistry |
| murcko_scaffold | string |  | Bemis–Murcko scaffold of the molecule |
| generic_scaffold | string |  | Generic (atom-type-agnostic) version of the scaffold |


### Source and Deployment
- **Source:** `Local`
- **Source Type:** `Internal`
- **DockerHub**: [https://hub.docker.com/r/ersiliaos/eos4k4f](https://hub.docker.com/r/ersiliaos/eos4k4f)
- **Docker Architecture:** `AMD64`, `ARM64`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos4k4f.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos4k4f.zip)

### Resource Consumption
- **Model Size (Mb):** `1`
- **Environment Size (Mb):** `784`
- **Image Size (Mb):** `741.53`

**Computational Performance (seconds):**
- 10 inputs: `27.77`
- 100 inputs: `18.32`
- 10000 inputs: `82.05`

### References
- **Source Code**: [https://github.com/chembl/ChEMBL_Structure_Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline)
- **Publication**: [https://doi.org/10.1186/s13321-020-00456-1](https://doi.org/10.1186/s13321-020-00456-1)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2025`
- **Ersilia Contributor:** [arnaucoma24](https://github.com/arnaucoma24)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [MIT](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos4k4f
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos4k4f
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
