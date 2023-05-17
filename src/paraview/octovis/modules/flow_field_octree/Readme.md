# Flow Field (Octree)

Specific version of the [Flow Field](../../../flow/modules/flow_field/Readme.md) (see [Flow plugin](../../../../../readme.md#flow-plugin)) to handle simulation data from a stellar merger.

## Input

The following input ports are available:

| Input  | Description                                                  | Type      | Remark |
| ------ | ------------------------------------------------------------ | --------- | ------ |
| Octree | Points (representing cells) of the octree, especially containing IDs to identify their position in the tree. | Point Set |        |
| Seed   | Input seed for stream-, streak- and path lines.              | Poly Data |        |
| Stars  | Identified star positions and respective properties.         | Point Set |        |

### Octree

The input points representing the octree contain the following data fields:

| Data field          | Description                                                  | Data   | Type        | Remark   |
| ------------------- | ------------------------------------------------------------ | ------ | ----------- | -------- |
| Point ID            | IDs indicating the cell within the octree.                   | Scalar | Point-based |          |
| Velocity            | Velocity.                                                    | Vector | Point-based |          |
| Minimum X           | Left domain boundary.                                        | Scalar | Point-based |          |
| Maximum X           | Right domain boundary.                                       | Scalar | Point-based |          |
| Minimum Y           | Lower domain boundary.                                       | Scalar | Point-based |          |
| Maximum Y           | Upper domain boundary.                                       | Scalar | Point-based |          |
| Minimum Z           | Back domain boundary.                                        | Scalar | Point-based |          |
| Maximum Z           | Front domain boundary.                                       | Scalar | Point-based |          |
| Time                | Time step of the data.                                       | Scalar | Point-based |          |
| Angular frequency   | Angular frequency representing the grid rotation.            | Scalar | Point-based |          |
| Star Classification | Classification of each cell, indicating its corresponding star. | Scalar | Point-based | Optional |
| Arrays              | Interpolate at particle positions and store values along stream- or pathlines. |        | Point-based | Optional |

### Stars

The input points representing the centers of the stars contain the following data fields:

| Data field    | Description                                                  | Data   | Type        | Remark |
| ------------- | ------------------------------------------------------------ | ------ | ----------- | ------ |
| Star Velocity | Velocity of the stars at the respective center of mass.      | Vector | Point-based |        |
| Star Rotation | Rotation of the stars around their center of mass, given as angular velocity. | Vector | Point-based |        |

## Parameters

The following parameters can be set by the user:

| Parameter                    | Description                                                  | Type    | Accepted values                                              | Default value          |
| ---------------------------- | ------------------------------------------------------------ | ------- | ------------------------------------------------------------ | ---------------------- |
| Method                       | Type of integral lines to generate.                          | Enum    | Streamlines, Streaklines, Pathlines                          | Streamlines            |
| Direction                    | Forward or reverse time integration.                         | Enum    | Forward, Reverse                                             | Forward                |
| Number of lines              | Number of advection steps.                                   | Integer | &gt;0                                                        | 20                     |
| Force fixed time step        | Force a fixed time step for streak- and path line computation. | Boolean | True, False                                                  | False                  |
| Fixed time step              | Time step size if fixed and for streamline computation.      | Float   | &gt;0                                                        | 1.0                    |
| Locality method              | Definition of locality, i.e., how to define the rotation of the stars. | Enum    | None, Rotation around origin, Star velocity, Star rotation, Rigid body | Rotation around origin |
| Force fixed frequency        | Override the frequency given in the dataset.                 | Boolean | True, False                                                  | False                  |
| Rotational frequency (omega) | Angular frequency in case of override.                       | Float   | &gt;0                                                        | 1.0                    |

## Output

The result of the computation are either stream-, streak or pathlines with the following data fields:

| Data field        | Description                                                  | Data   | Type        | Remark |
| ----------------- | ------------------------------------------------------------ | ------ | ----------- | ------ |
| ID (advection)    | Time information, storing the advection step.                | Scalar | Point-based |        |
| ID (distribution) | Time information, storing the advection step normalized for each integral line. | Scalar | Point-based |        |
| ID (advection)    | Time information, storing the advection step.                | Scalar | Cell-based  |        |
| ID (distribution) | Time information, storing the advection step normalized for each integral line. | Scalar | Cell-based  |        |

