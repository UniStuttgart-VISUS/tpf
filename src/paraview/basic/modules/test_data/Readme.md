# Test Data

Create simple 3D test geometry, represented as volume of fluid field and velocity field.

## Parameters

| Parameter                 | Description                                                                           | Type          | Accepted values   | Default value |
|---------------------------|---------------------------------------------------------------------------------------|---------------|-------------------|---------------|
| Geometry                  | Choose the output geometry as either *Sphere* or *Torus*.                             | Enumeration   | Sphere, Torus     | Sphere        |
| Number of cells           | Number of cells in each spatial direction for a cubic domain.                         | Integer       | > 0               | 50            |
| Rotation                  | Add velocity for the rotation around the objects's axis.                              | Boolean       | True, False       | False         |
| Rotation magnitude        | Velocity magnitude for the rotation.                                                  | Float         | > 0               | 1             |
| Translation               | Add velocity for the translation of the object.                                       | Boolean       | True, False       | False         |
| Translation magnitude     | Velocity magnitude for the translation.                                               | Float         | > 0               | 1             |
| Expansion                 | Add velocity for the expansion of the object.                                         | Boolean       | True, False       | False         |
| Expansion magnitude       | Velocity magnitude for the expansion.                                                 | Float         | > 0               | 1             |
| Torus rotation            | Add velocity for the rotation around the torus' center line.                          | Boolean       | True, False       | False         |
| Torus rotation magnitude  | Velocity magnitude for the torus' internal rotation.                                  | Float         | > 0               | 1             |

## Output

The output is a uniform grid, stored as rectilinear grid, comprising the following data fields:

| Data field                | Description                                                                           | Remark    |
|---------------------------|---------------------------------------------------------------------------------------|-----------|
| Volume of fluid           | The volume of fluid field representing the geometric object.                          |           |
| Velocity                  | Combined velocity of the selected velocity parts, such as rotation and translation.   |           |
