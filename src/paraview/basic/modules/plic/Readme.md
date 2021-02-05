# PLIC

Use piecewise linear interface calculation (PLIC) for the reconstruction of the fluid interface from volume of fluid (VOF) fields.

## Input

The following input ports are available:

| Input | Description                              | Type             | Remark |
| ----- | ---------------------------------------- | ---------------- | ------ |
| Grid  | Grid containing a volume of fluid field. | Rectilinear Grid |        |

### Grid

The grid contains the following data fields:

| Data field            | Description                                                  | Data   | Type       | Remark                                                       |
| --------------------- | ------------------------------------------------------------ | ------ | ---------- | ------------------------------------------------------------ |
| Volume of fluid field | A volume of fluid field, whose entries are in the range [0, 1]. | Scalar | Cell-based |                                                              |
| Interface gradients   | Gradients in interface cells.                                | Vector | Cell-based | Optional. If not given, the gradients are computed by the filter itself. |

## Parameters

The following parameters can be set by the user:

| Parameter                | Description                                                  | Type    | Accepted values | Default value |
| ------------------------ | ------------------------------------------------------------ | ------- | --------------- | ------------- |
| Number of iterations     | Number of iterations for the iterative approximation of the interface. | Integer | \>0             | 15            |
| Perturbation of vertices | Perturbation on the input vertices to prevent degenerate cases, which often occur in the initial time steps of the simulation. | Double  | &#8805;0        | 1e-5          |

## Output

The output is a polygonal mesh, for which the following data fields are provided:

| Data field | Description                                                  | Data   | Type       | Remark                                                       |
| ---------- | ------------------------------------------------------------ | ------ | ---------- | ------------------------------------------------------------ |
| Error      | Error, defined as difference between input volume of fluid and volume enclosed by the reconstructed interface. | Scalar | Cell-based | This field is only available when setting the TPF debug flag. |

