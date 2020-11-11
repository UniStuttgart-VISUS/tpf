# Correct VOF

Correct the volume of fluid (VOF) field by removing isolated non-zero cells.

## Input

| Input                     | Description                                                                           | Type              | Remark    |
|---------------------------|---------------------------------------------------------------------------------------|-------------------|-----------|
| Grid                      | Grid containing a volume of fluid field.                                              | Rectilinear Grid  |           |

### Grid

| Data field                | Description                                                                           | Remark    |
|---------------------------|---------------------------------------------------------------------------------------|-----------|
| Volume of fluid field     | A volume of fluid field, whose entries are in the range [0, 1].                       |           |

## Output

| Data field                | Description                                                                           | Remark    |
|---------------------------|---------------------------------------------------------------------------------------|-----------|
| Volume of fluid field     | The corrected volume of fluid field.                                                  |           |

The grid itself is not modified.
