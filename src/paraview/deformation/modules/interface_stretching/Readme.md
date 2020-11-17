# Interface Stretching

Compute interface stretching by means of a metric tensor [1,2,3,4], where the stretching value expresses change relative to one. A value larger than one indicates stretching, and a value smaller than one indicates contraction in the respective direction.

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

| Data field                     | Description                                                  | Data   | Type       | Remark |
| ------------------------------ | ------------------------------------------------------------ | ------ | ---------- | ------ |
| Stretching (area)              | Stretching expressed per area, where values >1 indicate stretching, and values <1 indicate contraction. | Scalar | Cell-based |        |
| Stretching direction (minimum) | Direction corresponding to minimum stretching.               | Vector | Cell-based |        |
| Stretching direction (maximum) | Direction corresponding to maximum stretching.               | Vector | Cell-based |        |
| Stretching direction (largest) | Direction corresponding to strongest stretching or bending.  | Vector | Cell-based |        |

The grid itself is not modified.

---

[1] Alexander Straub. Visualization of Interface Instabilities in Two-Phase Flow. *University of Stuttgart,* 2016.

[2] Alexander Straub, Grzegorz K. Karch, Sebastian Boblest, Jonas Kaufmann, Filip Sadlo, Bernhard Weigand, and Thomas Ertl. Visual Analysis of Interface Deformation in Multiphase Flow. *Proceedings of the DIPSI Workshop 2018,* *Università degli studi di Bergamo,* 45–47, 2018.

[3] Alexander Straub, Moritz Heinemann, and Thomas Ertl. Visualization and Visual Analysis for Multiphase Flow. *Proceedings of the DIPSI Workshop 2019,* *Università degli studi di Bergamo,* 25–27, 2019.

[4] Alexander Straub, and Thomas Ertl. Visualization Techniques for Droplet Interfaces and Multiphase Flow. *Droplet Interactions and Spray Processes,* *Springer International Publishing,* 121: 203–214, 2020.