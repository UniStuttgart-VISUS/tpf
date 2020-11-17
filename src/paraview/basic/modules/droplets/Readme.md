# Droplets

Extract droplets, defined as area of non-zero values in the volume of fluid (VOF) field. Additionally, properties such as volume and average velocity are calculated.

## Input

The following input ports are available:

| Input | Description                                                  | Type             | Remark |
| ----- | ------------------------------------------------------------ | ---------------- | ------ |
| Grid  | Grid containing fields representing fluids in a multiphase setting. | Rectilinear Grid |        |

### Grid

The grid contains the following data fields:

| Data field            | Description                                                  | Data   | Type       | Remark   |
| --------------------- | ------------------------------------------------------------ | ------ | ---------- | -------- |
| Volume of fluid field | A volume of fluid field, whose entries are in the range [0, 1]. | Scalar | Cell-based |          |
| Fluid positions       | Representative positions in cells containing the fluid phase. | Vector | Cell-based |          |
| Velocities            | Velocity field describing the fluid flow.                    | Vector | Cell-based | Optional |

## Parameters

The following parameters can be set by the user:

| Parameter                       | Description                                                  | Type        | Accepted values                           | Default value |
| ------------------------------- | ------------------------------------------------------------ | ----------- | ----------------------------------------- | ------------- |
| Calculate translation           | Calculate the average velocity for each droplet.             | Boolean     |                                           | False         |
| Calculate rotation              | Calculate the rotation of the droplet, assuming rigid-body dynamics. | Boolean     |                                           | False         |
| Calculate energy                | Calculate the energy related to droplet movement.            | Boolean     |                                           | False         |
| Calculate inertia               | Calculate the inertia tensor for the droplet                 | Boolean     |                                           | False         |
| Inertia scaling                 | Scale the inertia vectors (see output below) for visualization output by *normalizing*, *eigenvalue*, or *logarithm of the eigenvalue*. | Enumeration | Normalized, Eigenvalue, Log of Eigenvalue | Normalized    |
| Method for rotation computation | Method for rotation computation: *fluid mechanics*, *velocities*, or *PCA*. | Enumeration | Mechanics, Velocities, PCA                | Mechanics     |

## Output

The following output ports are available:

| Input     | Description                                                  | Type             | Remark |
| --------- | ------------------------------------------------------------ | ---------------- | ------ |
| Grid      | A copy of the input grid.                                    | Rectilinear Grid |        |
| Positions | A set of points representing the center of mass of the droplets. | Polydata         |        |

### Grid

The following data fields are appended to the input grid:

| Data field | Description                                              | Data   | Type       | Remark |
| ---------- | -------------------------------------------------------- | ------ | ---------- | ------ |
| Droplet ID | Assigned ID of the droplet, to match grid and positions. | Scalar | Cell-based |        |
| Volume     | Volume of the droplet.                                   | Scalar | Cell-based |        |

The grid itself is not modified.

### Positions

The following data fields are appended to the positions:

| Data field         | Description                                                  | Data   | Type        | Remark                                                       |
| ------------------ | ------------------------------------------------------------ | ------ | ----------- | ------------------------------------------------------------ |
| Droplet ID         | Assigned ID of the droplet, to match grid and positions.     | Scalar | Point-based |                                                              |
| Volume             | Volume of the droplet.                                       | Scalar | Point-based |                                                              |
| Radius             | Representative radius of the droplet for visualization.      | Scalar | Point-based |                                                              |
| Translation        | Average velocity of the droplet.                             | Vector | Point-based | Requires velocities.<br />Translation computation must be enabled. |
| Translation energy | Kinetic energy of the average droplet velocity.              | Scalar | Point-based | Requires velocities.<br />Translation, and energy computation must be enabled. |
| Rotation           | Scaled rotation axis, defining the rotation of the droplet.  | Vector | Point-based | Requires velocities.<br />Rotation computation must be enabled. |
| Rotation energy    | Kinetic energy of the droplet rotation.                      | Scalar | Point-based | Requires velocities.<br />Rotation, and energy computation must be enabled. |
| Local energy       | Internal energy of the droplet, which is defined as the total energy subtracting translation and rotation energy. | Scalar | Point-based | Requires velocities.<br />Translation, rotation, and energy computation must be enabled. |
| Energy             | Total energy of the droplet, defined as kinetic energy of the input vector field. | Scalar | Point-based | Requires velocities.<br />Energy computation must be enabled. |
| Inertia            | Inertia tensor for the droplet.                              | Matrix | Point-based | Requires velocities.<br />Inertia computation must be enabled. |
| Inertia #1         | First eigenvector of the inertia tensor, scaled by its corresponding eigenvalue according to the scaling method. | Vector | Point-based | Requires velocities.<br />Inertia computation, and detailed output must be enabled. |
| Inertia #2         | Second eigenvector of the inertia tensor, scaled by its corresponding eigenvalue according to the scaling method. | Vector | Point-based | Requires velocities.<br />Inertia computation, and detailed output must be enabled. |
| Inertia #3         | Third eigenvector of the inertia tensor, scaled by its corresponding eigenvalue according to the scaling method. | Vector | Point-based | Requires velocities.<br />Inertia computation, and detailed output must be enabled. |