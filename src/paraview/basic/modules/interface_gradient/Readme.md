# Interface Gradient

Calculate the gradient in interface cells, i.e., in cells where the volume of fluid (VOF) value is larger than 0 and smaller than 1.

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

## Output

The following data fields are appended to the input grid:

| Data field         | Description                                                  | Data   | Type       | Remark |
| ------------------ | ------------------------------------------------------------ | ------ | ---------- | ------ |
| Interface gradient | Within interface cells, it contains the computed gradient, elsewhere it is set to a zero-vector. | Vector | Cell-based |        |

The grid itself is not modified.

