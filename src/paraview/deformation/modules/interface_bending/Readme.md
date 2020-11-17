# Interface Bending

Compute interface bending as change in curvature [1,2,3], where the bending value expresses this change relative to zero. A negative value indicates an increase in convexity (or decrease in concavity), and a positive value an increase in concavity (or decrease in convexity) in the respective direction.

## Input

The following input ports are available:

| Input | Description                                                  | Type             | Remark |
| ----- | ------------------------------------------------------------ | ---------------- | ------ |
| Grid  | Grid containing fields representing fluids in a multiphase setting. | Rectilinear Grid |        |

### Grid

The grid contains the following data fields:

| Data field            | Description                                                  | Data   | Type       | Remark |
| --------------------- | ------------------------------------------------------------ | ------ | ---------- | ------ |
| Volume of fluid field | A volume of fluid field, whose entries are in the range [0, 1]. | Scalar | Cell-based |        |
| Interface gradients   | Gradient at interface cells.                                 | Vector | Cell-based |        |
| Interface positions   | Field containing the interface barycenter of each interface cell. | Vector | Cell-based |        |
| Velocities            | Velocity field describing the fluid flow.                    | Vector | Cell-based |        |

## Output

The following data fields are appended to the input grid:

| Data field                           | Description                                                  | Data   | Type       | Remark |
| ------------------------------------ | ------------------------------------------------------------ | ------ | ---------- | ------ |
| Bending (minimum)                    | Minimum bending, where bending values >0 indicate increase in concavity, and values <0 indicate increase in convexity. | Scalar | Cell-based |        |
| Bending (maximum)                    | Maximum bending, where bending values >0 indicate increase in concavity, and values <0 indicate increase in convexity. | Scalar | Cell-based |        |
| Bending (absolute maximum)           | Strongest bending, where bending values >0 indicate increase in concavity, and values <0 indicate increase in convexity. | Scalar | Cell-based |        |
| Bending direction (minimum)          | Direction corresponding to minimum bending.                  | Vector | Cell-based |        |
| Bending direction (maximum)          | Direction corresponding to maximum bending.                  | Vector | Cell-based |        |
| Bending direction (absolute maximum) | Direction corresponding to strongest bending.                | Vector | Cell-based |        |

The grid itself is not modified.

---

[1] Alexander Straub, Grzegorz K. Karch, Sebastian Boblest, Jonas Kaufmann, Filip Sadlo, Bernhard Weigand, and Thomas Ertl. Visual Analysis of Interface Deformation in Multiphase Flow. *Proceedings of the DIPSI Workshop 2018,* *Università degli studi di Bergamo,* 45–47, 2018.

[2] Alexander Straub, Moritz Heinemann, and Thomas Ertl. Visualization and Visual Analysis for Multiphase Flow. *Proceedings of the DIPSI Workshop 2019,* *Università degli studi di Bergamo,* 25–27, 2019.

[3] Alexander Straub, and Thomas Ertl. Visualization Techniques for Droplet Interfaces and Multiphase Flow. *Droplet Interactions and Spray Processes,* *Springer International Publishing,* 121: 203–214, 2020.