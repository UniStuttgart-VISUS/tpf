# Binary Star

Classification of binary stars in a stellar merger, defined by physical properties.

## Input

The following input ports are available:

| Input  | Description                                                  | Type      | Remark |
| ------ | ------------------------------------------------------------ | --------- | ------ |
| Octree | Points (representing cells) of the octree, especially containing IDs to identify their position in the tree. | Point Set |        |

### Octree

The input points representing the octree contain the following data fields:

| Data field         | Description                                                  | Data   | Type        | Remark |
| ------------------ | ------------------------------------------------------------ | ------ | ----------- | ------ |
| Density            | Local density of the star.                                   | Scalar | Point-based |        |
| Density (Accretor) | Density of the phase belonging to the accretor.              | Scalar | Point-based |        |
| Density (Donor)    | Density of the phase belonging to the donor.                 | Scalar | Point-based |        |
| Density (Other)    | Density of the phase belonging to neither accretor nor donor. | Scalar | Point-based |        |
| Velocity           | Velocity.                                                    | Vector | Point-based |        |
| Gravitation        | Gravitation towards the respective center of the star.       | Vector | Point-based |        |
| Cell size          | Sizes of the octree cells.                                   | Scalar | Point-based |        |
| Minimum X          | Left domain boundary.                                        | Scalar | Point-based |        |
| Maximum X          | Right domain boundary.                                       | Scalar | Point-based |        |
| Minimum Y          | Lower domain boundary.                                       | Scalar | Point-based |        |
| Maximum Y          | Upper domain boundary.                                       | Scalar | Point-based |        |
| Minimum Z          | Back domain boundary.                                        | Scalar | Point-based |        |
| Maximum Z          | Front domain boundary.                                       | Scalar | Point-based |        |

## Parameters

The following parameters can be set by the user:

| Parameter            | Description                                                  | Type    | Accepted values | Default value |
| -------------------- | ------------------------------------------------------------ | ------- | --------------- | ------------- |
| Number of iterations | Number of iteration steps for the algorithm.                 | Integer | &gt;0           | 5             |
| Density cutoff       | Cells with a density below this cutoff are discarded and not considered to be part of the stars. | Float   | &geq;0          | 1.0           |

## Output

The following output ports are available:

| Input  | Description                                                  | Type      | Remark |
| :----- | :----------------------------------------------------------- | :-------- | :----- |
| Octree | A copy of the input octree, enriched with classification and information for each octree cell. | Point Set |        |
| Stars  | Points representing the two stars, together with some properties. | Poly Data |        |

### Octree

The following data fields are appended to the input octree:

| Data field        | Description                                                  | Data   | Type        | Remark |
| ----------------- | ------------------------------------------------------------ | ------ | ----------- | ------ |
| Classification    | Classified as belonging to the accretor (1), the donor (2), or the surrounding phase (0). | Scalar | Point-based |        |
| Orbital frequency | Angular frequency of the rotation of the stars around a common center (of mass). | Scalar | Global      |        |

The points of the octree and the input data fields are not modified.

### Stars

These are the points representing the centers of the stars, together with these data fields:

| Data field        | Description                                                  | Data   | Type        | Remark |
| ----------------- | ------------------------------------------------------------ | ------ | ----------- | ------ |
| Velocity          | Velocity at the center of mass of the respective star.       | Vector | Point-based |        |
| Mass              | Mass of the respective star.                                 | Scalar | Point-based |        |
| Angular velocity  | Angular velocity of the respective star.                     | Vector | Point-based |        |
| Orbital frequency | Angular frequency of the rotation of the stars around a common center (of mass). | Scalar | Point-based |        |
| Roche lobe radius | Roche lobe radius.                                           | Scalar | Point-based |        |
