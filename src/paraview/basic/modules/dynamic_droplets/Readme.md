# Dynamic droplets

Create static glyphs representing the rotation and translation of the time-dependent input droplets.

## Input

The following input ports are available:

| Input | Description                                                           | Type              | Remark |
| ----- | --------------------------------------------------------------------- | ----------------- | ------ |
| Input | Extracted droplet information, such as translation and rotation.      | Poly Data         |        |
| Grid  | Grid containing fields representing fluids in a multiphase setting.   | Rectilinear Grid  |        |

### Input

The input points representing the droplet positions contain the following additional data fields:

| Data field    | Description                                       | Data   | Type         | Remark   |
| ------------- | ------------------------------------------------- | ------ | ------------ | -------- |
| Translation   | Translation of the droplet.                       | Vector | Point-based  |          |
| Rotation      | Rotation axis of the droplet.                     | Vector | Point-based  |          |
| Radius        | Approximate representative radius of the droplet. | Scalar | Point-based  |          |

### Grid

The grid contains the following data fields:

| Data field    | Description                       | Data   | Type       | Remark   |
| ------------- | --------------------------------- | ------ | ---------- | -------- |
| Droplet IDs   | IDs of the extracted droplets.    | Scalar | Cell-based |          |

## Parameters

The following parameters can be set by the user:

| Parameter                         | Description                                                  | Type       | Accepted values                           | Default value |
| --------------------------------- | ------------------------------------------------------------ | ---------- | ----------------------------------------- | ------------- |
| Number of time steps              | Number of time steps for which information is represented.   | Integer    | \>0                                       | 10            |
| Compute                           | Switch to allow not to compute information directly.         | Boolean    |                                           | False         |
| Static frame of reference         | Use the initial droplet information only.                    | Boolean    |                                           | False         |
| Ribbon size                       | Width scalar for the created ribbons.                        | Float      | \>0                                       | 0.025         |
| Fix size of rotation axis         | Fix size of the rotation axis to the circumsphere.           | Boolean    |                                           | True          |
| Rotation axis scale               | Scale rotation axis if size is not fixed.                    | Float      | \>0                                       | 1.0           |
| Show rotation at the origin       | Show rotation axes at the original droplet's origin.         | Boolean    |                                           | True          |

## Output

The output is a multiblock dataset consisting of several blocks:

| Block name            | Description                                                       | Type           | Remark            |
| --------------------- | ----------------------------------------------------------------- | -------------- | ----------------- |
| Droplet tracks        | Line strip showing the movement of the centers of mass.           | Poly Data      |                   |
| Translation summary   | Point on the circumsphere in the direction of the translation.    | Poly Data      |                   |
| Translation paths     | Line strip beginning on the circumsphere, showing translation.    | Poly Data      |                   |
| Rotation axes         | Lines representing a rotation axis per time step.                 | Poly Data      |                   |
| Rotation ribbons      | Ribbons from connecting translated rotation axes.                 | Poly Data      |                   |
| Rotation paths        | Bent arrow glyph along the surface of the circumsphere.           | Poly Data      |                   |
| Coordinate axes       | Ribbons from connecting coordinate system axes (experimental).    | Poly Data      |                   |

The following data fields are attached to the blocks (see remark for availability):

| Data field            | Description                                                               | Data   | Type                 | Remark                                                                            |
| --------------------- | ------------------------------------------------------------------------- | ------ | -------------------- | --------------------------------------------------------------------------------- |
| Topology information  | Droplet topology information: 0 - unchanged, 1 - breakup, 2 - collision.  | Scalar | Cell-based           | Available for: Droplet tracks                                                     |
| Displacement          | Displacement vector indicating total droplet translation.                 | Vector | Point-based          | Available for: Translation summary                                                |
| Time ID               | Relative integer time step.                                               | Scalar | Point- or cell-based | Available for: Translation paths, rotation axes, rotation ribbons, rotation paths |
| Time information      | Time: [0, 1] on x-axis, [2, 3] on y-axis and [4, 5] on z-axis.            | Scalar | Point-based          | Available for: Coordinate axes                                                    |
