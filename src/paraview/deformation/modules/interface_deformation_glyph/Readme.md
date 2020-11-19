# Interface Deformation Glyph

Provides a glyph specialized for visualizing interface deformation.

## Input

The following input ports are available:

| Input | Description                                                  | Type             | Remark |
| ----- | ------------------------------------------------------------ | ---------------- | ------ |
| Grid  | Grid containing fields representing fluids in a multiphase setting. | Rectilinear Grid |        |

### Grid

The grid contains the following data fields:

| Data field                     | Description                                                  | Data   | Type       | Remark |
| ------------------------------ | ------------------------------------------------------------ | ------ | ---------- | ------ |
| Interface positions            | Field containing the interface barycenter of each interface cell. | Vector | Cell-based |        |
| Velocities                     | Velocity field describing the fluid flow.                    | Vector | Cell-based | Optional |
| Stretching direction (minimum) | Direction corresponding to minimum stretching.               | Vector | Cell-based | Optional |
| Stretching direction (maximum) | Direction corresponding to maximum stretching.               | Vector | Cell-based | Optional |
| Bending (minimum)                    | Minimum bending, where bending values >0 indicate increase in concavity, and values <0 indicate increase in convexity. | Scalar | Cell-based | Optional |
| Bending (maximum)                    | Maximum bending, where bending values >0 indicate increase in concavity, and values <0 indicate increase in convexity. | Scalar | Cell-based | Optional |
| Bending direction (minimum)          | Direction corresponding to minimum bending.                  | Vector | Cell-based | Optional |
| Bending direction (maximum)          | Direction corresponding to maximum bending.                  | Vector | Cell-based | Optional |

The data fields can be considered as optional, however the following constellations are required for different features:

| Feature          | Required data fields                                         |
| ---------------- | ------------------------------------------------------------ |
| Velocity glyph   | Velocities                                                   |
| Stretching glyph | Stretching direction (minimum and maximum)                   |
| Bending glyph    | Bending (minimum and maximum), Bending direction (minimum and maximum) |

## Output

The output is polygonal mesh, represented by a polydata object. It contains the following data fields useful for coloring the glyph:

| Data field         | Description                                            | Data   | Type         | Remark |
| ------------------ | ------------------------------------------------------ | ------ | ------------ | ------ |
| Velocity magnitude | Magnitude of the velocity.                             | Scalar | Object-based |        |
| Stretching (area)  | Combined stretching for minimum and maximum (product). | Scalar | Object-based |        |
| Bending            | Combined bending for minimum and maximum (sum).        | Scalar | Object-based |        |

