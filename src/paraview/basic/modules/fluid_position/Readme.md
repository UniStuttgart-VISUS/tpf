# Fluid Position

Calculate a representative position for the fluid in a cell. There are three positions that can be computed: cell center, fluid center of mass, and interface barycenter.

## Input

The following input ports are available:

| Input | Description                                                  | Type             | Remark |
| ----- | ------------------------------------------------------------ | ---------------- | ------ |
| Grid  | Grid containing fields representing fluids in a multiphase setting. | Rectilinear Grid |        |

### Grid

The grid contains the following data fields:

| Data field            | Description                                                  | Data   | Type       | Remark |
| --------------------- | ------------------------------------------------------------ | ------ | ---------- | ------ |
| Volume of fluid field | A volume of fluid field, whose entries are in the range [0, 1]. | Scalar | Cell-based |        |
| Interface gradients   | Gradients in interface cells.                                | Vector | Cell-based |        |

## Parameters

The following parameters can be set by the user:

| Parameter     | Description                                                  | Type        | Accepted values                      | Default value |
| ------------- | ------------------------------------------------------------ | ----------- | ------------------------------------ | ------------- |
| Position type | Choose the position type as either *cell center*, *fluid center of mass*, or *interface barycenter*. | Enumeration | Cell center, Fluid center, Interface | Cell center   |

## Output

The following output ports are available:

| Input     | Description                                            | Type             | Remark |
| --------- | ------------------------------------------------------ | ---------------- | ------ |
| Grid      | A copy of the input grid.                              | Rectilinear Grid |        |
| Positions | A set of points representing the calculated positions. | Polydata         |        |

### Grid

The following data fields are appended to the input grid:

| Data field     | Description                                             | Data   | Type       | Remark |
| -------------- | ------------------------------------------------------- | ------ | ---------- | ------ |
| Fluid position | The representative position of the fluid within a cell. | Vector | Cell-based |        |

The grid itself is not modified.

### Positions

The points representing the positions do not have any attributes attached.