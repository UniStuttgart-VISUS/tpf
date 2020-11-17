# Correct VOF

Correct the volume of fluid (VOF) field by removing isolated non-zero cells.

## Input

The following input ports are available:

| Input                     | Description                                                                           | Type              | Remark    |
|---------------------------|---------------------------------------------------------------------------------------|-------------------|-----------|
| Grid                      | Grid containing a volume of fluid field.                                              | Rectilinear Grid  |           |

### Grid

The grid contains the following data fields:

| Data field            | Description                                                  | Data   | Type       | Remark |
| --------------------- | ------------------------------------------------------------ | ------ | ---------- | ------ |
| Volume of fluid field | A volume of fluid field, whose entries are in the range [0, 1]. | Scalar | Cell-based |        |

## Output

The following data fields are appended to the input grid:

| Data field            | Description                          | Data   | Type       | Remark |
| --------------------- | ------------------------------------ | ------ | ---------- | ------ |
| Volume of fluid field | The corrected volume of fluid field. | Scalar | Cell-based |        |

The grid itself is not modified.

