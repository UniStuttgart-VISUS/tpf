# Compute Temperature

Compute the temperature from the energies and densities.

## Input

The following input ports are available:

| Input | Description                               | Type     | Remark |
| ----- | ----------------------------------------- | -------- | ------ |
| Data  | Points containing energies and densities. | Data Set |        |

### Data

The input points must contain the following data fields:

| Data field       | Description                                | Data   | Type          | Remark |
| ---------------- | ------------------------------------------ | ------ | ------------- | ------ |
| Density          | Local density of the star.                 | Scalar | Point-based   |        |
| Density First    | Density of the first specie.               | Scalar | Point-based   |        |
| Density Second   | Density of the second specie.              | Scalar | Point-based   |        |
| Density Third    | Density of the third specie.               | Scalar | Point-based   |        |
| Density Fourth   | Density of the fourth specie.              | Scalar | Point-based   |        |
| Density Fifth    | Density of the fifth specie.               | Scalar | Point-based   |        |
| Total Energy     | Total energy (E).                          | Scalar | Point-based   |        |
| Momentum         | Inertial momentum.                         | Vector | Point-based   |        |
| Tracer           | Entropy tracer (tau).                      | Scalar | Point-based   |        |
| Molecular Masses | The molecular masses for all five species. | Vector | Dataset-based |        |
| Atomic Numbers   | The atomic numbers for all five species.   | Vector | Dataset-based |        |

## Output

Copy of the input with added temperature field T.

