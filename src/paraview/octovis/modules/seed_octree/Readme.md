# Seed Octree

Use the octree data structure to create a seed.

## Input

The following input ports are available:

| Input  | Description                                                  | Type      | Remark |
| ------ | ------------------------------------------------------------ | --------- | ------ |
| Octree | Points (representing cells) of the octree, especially containing IDs to identify their position in the tree. | Point Set |        |

### Octree

The input points representing the octree contain the following data fields:

| Data field | Description                                | Data   | Type        | Remark   |
| ---------- | ------------------------------------------ | ------ | ----------- | -------- |
| Point ID   | IDs indicating the cell within the octree. | Scalar | Point-based |          |
| Density    | Local density of the star.                 | Scalar | Point-based | Optional |
| Minimum X  | Left domain boundary.                      | Scalar | Point-based |          |
| Maximum X  | Right domain boundary.                     | Scalar | Point-based |          |
| Minimum Y  | Lower domain boundary.                     | Scalar | Point-based |          |
| Maximum Y  | Upper domain boundary.                     | Scalar | Point-based |          |
| Minimum Z  | Back domain boundary.                      | Scalar | Point-based |          |
| Maximum Z  | Front domain boundary.                     | Scalar | Point-based |          |

## Parameters

The following parameters can be set by the user:

| Parameter   | Description                                                  | Type    | Accepted values           | Default value |
| ----------- | ------------------------------------------------------------ | ------- | ------------------------- | ------------- |
| Seed method | Seeding strategy.                                            | Enum    | Octree leaves, Isosurface | Octree leaves |
| Seed offset | Start seeding only in the n-th cell.                         | Integer | &geq;0                    | 0             |
| Seed size   | Only seed in a limited number of cells.                      | Integer | &geq;0 and -1 (ignore)    | -1            |
| Isovalue    | Isovalue for density-based seeding at the isosurface. *Requires density input to be available.* | Float   |                           | 0.5           |
| Predicate   | Only seed in cells that contain the isosurface if the predicate is fulfilled for the cell center (seed position). | Enum    | Both, Larger, Smaller     | Both          |

## Output

The output consists of the seed points.

