# FS3D Reader

Read temporal FS3D files (.lst) containing grid information, together with either a scalar field or a vector field.

## File specification

The information is spread into multiple files, where only the *.lst* file is given as input parameter:

| File (type)               | Description                                                                                   | Remark            |
|---------------------------|-----------------------------------------------------------------------------------------------|-------------------|
| .lst                      | This is an ASCII file containing a list of temporally subsequent input files, one per line.   |                   |
| netz_xyz.\[bin\|dat\]     | A binary file storing the spatial information of the rectilinear grid.                        |                   |
| .\[bin\|dat\]             | Binary files storing the data fields and respective physical information per time step.       |                   |

Please note that the *netz_xyz.\[bin\|dat\]* file has to be stored next to the input *.lst* file. Only files listed in the *.lst* file are loaded, their paths being relative to the *.lst* file itself.

### netz_xyz.\[bin\|dat\]

Every file contains the following header, with offset and length given in bytes:

| Offset    | Length    | Type      | Description                                                           |
|-----------|-----------|-----------|-----------------------------------------------------------------------|
| 0         | 80        | chars     | Unit of length as human-readable string.                              |
| 80        | 4         | int       | Resolution in x direction.                                            |
| 84        | 4         | int       | Resolution in y direction.                                            |
| 88        | 4         | int       | Resolution in z direction.                                            |

The data is stored as follows:

| Offset                                    | Length            | Type      | Description                                                           |
|-------------------------------------------|-------------------|-----------|-----------------------------------------------------------------------|
| 92                                        | x-resolution * 8  | double    | x components of the stored positions.                                 |
| 92 + x-resolution * 8                     | y-resolution * 8  | double    | y components of the stored positions.                                 |
| 92 + (x-resolution + y-resolution) * 8    | z-resolution * 8  | double    | z components of the stored positions.                                 |

### .\[bin\|dat\] files

Every file contains the following header, with offset and length given in bytes:

| Offset    | Length    | Type      | Description                                                                           |
|-----------|-----------|-----------|---------------------------------------------------------------------------------------|
| 0         | 80        | chars     | Name of the stored data field as human-readable string.                               |
| 80        | 80        | chars     | Unit of the stored data field as human-readable string.                               |
| 160       | 4         | int       | The ID of the time step as sequential number starting at 0 for the first time step.   |
| 164       | 8         | double    | The simulation time.                                                                  |
| 172       | 4         | int       | Resolution in x direction.                                                            |
| 176       | 4         | int       | Resolution in y direction.                                                            |
| 180       | 4         | int       | Resolution in z direction.                                                            |

The data is stored in z-major ordering for ___scalar values___ as follows:

| Offset                            | Length                                            | Type      | Description                   |
|-----------------------------------|---------------------------------------------------|-----------|-------------------------------|
| 188                               | x-resolution * y-resolution * z-resolution * 8    | double    | Scalar values.                |

The data is stored in z-major ordering for ___vectors___ as follows:

| Offset                                                | Length                                            | Type      | Description                   |
|-------------------------------------------------------|---------------------------------------------------|-----------|-------------------------------|
| 188                                                   | x-resolution * y-resolution * z-resolution * 8    | double    | x components of the vectors.  |
| 188 + x-resolution * y-resolution * z-resolution * 8  | x-resolution * y-resolution * z-resolution * 8    | double    | y components of the vectors.  |
| 188 + x-resolution * y-resolution * z-resolution * 16 | x-resolution * y-resolution * z-resolution * 8    | double    | z components of the vectors.  |

## Output

The output is a rectilinear grid with either one of the following data fields:

| Data field                | Description                       | Remark                                                                    |
|---------------------------|-----------------------------------|---------------------------------------------------------------------------|
| Scalar field              | A scalar field.                   | A scalar field is assumed when the field stores a single value per node.  |
| Vector field              | A three-dimensional vector field. | A vector field is assumed when the field stores three values per node.    |

The name of the respective output field is determined by the name in the input file.
