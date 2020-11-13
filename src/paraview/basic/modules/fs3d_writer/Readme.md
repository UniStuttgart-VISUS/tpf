# FS3D Writer

Save a rectilinear grid and attached data field to FS3D files.

## File specification

Please see the documentation of the [FS3D Reader](../fs3d_reader/Readme.md) for detailed specification of the output file format.

## Input

The following input ports are available:

| Input | Description                  | Type             | Remark |
| ----- | ---------------------------- | ---------------- | ------ |
| Grid  | Grid containing data fields. | Rectilinear Grid |        |

### Grid

The grid has to contain at least one data field, which is supposed to be saved to file.

## Parameters

The following parameters can be set by the user:

| Parameter    | Description                                     | Type            | Accepted values | Default value |
| ------------ | ----------------------------------------------- | --------------- | --------------- | ------------- |
| Output field | The data field name, which is written to file.  | Data field name |                 |               |
| File name    | Output file name.                               | Path            |                 |               |
| Name         | Name of the data field, describing its content. | String          |                 | Data          |
| Unit         | Unit of the stored data.                        | String          |                 | [-]           |
| Grid unit    | Length unit of the underlying grid.             | String          |                 | [cm]          |

## Output

The input grid is passed through without modification.

