# Flow Field

Stream-, streak- and pathline computation, with emphasis on per-droplet frames of reference.

## Input

The following input ports are available:

| Input    | Description                                                  | Type             | Remark |
| -------- | ------------------------------------------------------------ | ---------------- | ------ |
| Input    | Points (representing cells) of the octree, especially containing IDs to identify their position in the tree. | Rectilinear Grid |        |
| Seed     | Input seed for stream-, streak- and path lines.              | Poly Data        |        |
| Droplets | Extracted droplet positions with corresponding properties.   | Poly Data        |        |

### Input

The input rectilinear grid contains the following data fields:

| Data field            | Description                                                  | Data   | Type        | Remark   |
| --------------------- | ------------------------------------------------------------ | ------ | ----------- | -------- |
| Droplet ID            | ID to relate the cell to the corresponding droplet.          | Scalar | Cell-based  | Optional |
| Velocity field        | Velocity at the cell's center.                               | Vector | Cell-based  |          |
| Volume of fluid field | Volume of fluid field, giving the volume fraction for each cell. | Scalar | Cell-based  | Optional |
| Arrays                | Interpolate at particle positions and store values along stream- or pathlines. |        | Point-based | Optional |

### Droplets

The input points representing the centers of the droplets contain the following data fields:

| Data field          | Description                                                  | Data   | Type        | Remark |
| ------------------- | ------------------------------------------------------------ | ------ | ----------- | ------ |
| Droplet translation | Translation of the droplet at its center of mass.            | Vector | Point-based |        |
| Droplet rotation    | Rotation of the droplet around its axis of rotation, given as angular velocity. | Vector | Point-based |        |

## Parameters

The following parameters can be set by the user:

| Parameter            | Description                                                  | Type           | Accepted values                           | Default value   |
| -------------------- | ------------------------------------------------------------ | -------------- | ----------------------------------------- | --------------- |
| Seed in cells        | Seed in cells instead of using the input seed.               | Boolean        | True, False                               | True            |
| Cell type            | Choose in which cells seeds should be created.               | Enum           | All cells, Interface cells, Droplet cells | All cells       |
| Seed size per cell   | Number of seed points per cell.                              | Integer        | &gt;0                                     | 1               |
| Method               | Type of the integral curves to compute.                      | Enum           | Streamlines, Streaklines, Pathlines       | Streamlines     |
| Interpolatable       | Is the data interpolatable? (Not the case for VOF data.)     | Boolean        | True, False                               | False           |
| Integration          | Integration method. Availability is based on interpolatability of the input. | Enum           | Explicit Euler, Adams-Bashforth           | Adams-Bashforth |
| Number of advections | Number of advection steps for streamline computation.        | Integer        | &gt;0                                     | 20              |
| Time range           | Time range for pathline and streakline computation.          | [Float, Float] |                                           | [0.0, 1.0]      |
| Time step            | Time step size for streamline computation and if the data is interpolatable. | Float          |                                           | 0.1             |
| Time dependency      | Time dependency of the frame of reference, i.e., static vs. dynamic frame of reference. | Enum           | Dynamic, Static                           | Dynamic         |
| Keep translation     | Do not remove translational velocity part.                   | Boolean        | True, False                               | False           |
| Keep rotation        | Do not remove rotational velocity part.                      | Boolean        | True, False                               | False           |

## Output

The result of the computation are either stream-, streak or pathlines with the following data fields:

| Data field        | Description                                                  | Data   | Type        | Remark |
| ----------------- | ------------------------------------------------------------ | ------ | ----------- | ------ |
| ID (advection)    | Time information, storing the advection step.                | Scalar | Point-based |        |
| ID (distribution) | Time information, storing the advection step normalized for each integral line. | Scalar | Point-based |        |
| Time              | Time information, storing the time step.                     | Scalar | Point-based |        |
| Seed ID           | ID of the seed point, to discern individual integral lines.  | Scalar | Point-based |        |

