# Interface Deformation

Compute interface deformation by means of the derived quantities [interface stretching](../interface_stretching/Readme.md) [1,2,3,4] and [interface bending](../interface_bending/Readme.md) [2,3,4].

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
| Velocities            | Velocity field describing the fluid flow.                    | Vector | Cell-based | Optional |

## Parameters

The following parameters can be set by the user:

| Parameter       | Description                                                  | Type    | Accepted values | Default value |
| --------------- | ------------------------------------------------------------ | ------- | --------------- | ------------- |
| Surface tension | Calculate and use surface tension force instead of input velocity field. Must be set to *True* when omitting the velocity field input.<br />*Note that this is an experimental feature and is probably far from physically correct.* | Boolean |                 | False         |
| Coefficient     | Surface tension coefficient based on the involved species.   | Float   | \>0             | 72.75         |
| Density         | Density of the liquid phase.                                 | Float   | \>0             | 0.9982        |
| Time step       | Alternative time step for scaling the surface tension force vectors. If set to zero, use time step derived from the dataset. | Float   | &#8805;0        | 0.0           |

## Output

The following data fields are appended to the input grid:

| Data field                           | Description                                                  | Data   | Type       | Remark                                                       |
| ------------------------------------ | ------------------------------------------------------------ | ------ | ---------- | ------------------------------------------------------------ |
| Stretching (area)                    | Stretching expressed per area, where values >1 indicate stretching, and values <1 indicate contraction. | Scalar | Cell-based |                                                              |
| Stretching direction (minimum)       | Direction corresponding to minimum stretching.               | Vector | Cell-based |                                                              |
| Stretching direction (maximum)       | Direction corresponding to maximum stretching.               | Vector | Cell-based |                                                              |
| Stretching direction (largest)       | Direction corresponding to strongest stretching or bending.  | Vector | Cell-based |                                                              |
| Bending (minimum)                    | Minimum bending, where bending values >0 indicate increase in concavity, and values <0 indicate increase in convexity. | Scalar | Cell-based |                                                              |
| Bending (maximum)                    | Maximum bending, where bending values >0 indicate increase in concavity, and values <0 indicate increase in convexity. | Scalar | Cell-based |                                                              |
| Bending (absolute maximum)           | Strongest bending, where bending values >0 indicate increase in concavity, and values <0 indicate increase in convexity. | Scalar | Cell-based |                                                              |
| Bending direction (minimum)          | Direction corresponding to minimum bending.                  | Vector | Cell-based |                                                              |
| Bending direction (maximum)          | Direction corresponding to maximum bending.                  | Vector | Cell-based |                                                              |
| Bending direction (absolute maximum) | Direction corresponding to strongest bending.                | Vector | Cell-based |                                                              |
| Interface gradient                   | Within interface cells, it contains the computed gradient, elsewhere it is set to a zero-vector. | Vector | Cell-based | Only available when compiled with detailed output.           |
| Interface position                   | Field containing the interface barycenter of each interface cell. | Vector | Cell-based | Only available when compiled with detailed output.           |
| Interface curvature                  | The curvature at the interface barycenter.                   | Scalar | Cell-based | Only available when compiled with detailed output, and surface tension force computation turned on. |
| Surface tension force                | Approximate surface tension force in interface cells.        | Vector | Cell-based | Only available when compiled with detailed output, and surface tension force computation turned on. |

The grid itself is not modified.

---

[1] Alexander Straub. Visualization of Interface Instabilities in Two-Phase Flow. *University of Stuttgart,* 2016.

[2] Alexander Straub, Grzegorz K. Karch, Sebastian Boblest, Jonas Kaufmann, Filip Sadlo, Bernhard Weigand, and Thomas Ertl. Visual Analysis of Interface Deformation in Multiphase Flow. *Proceedings of the DIPSI Workshop 2018,* *Università degli studi di Bergamo,* 45–47, 2018.

[3] Alexander Straub, Moritz Heinemann, and Thomas Ertl. Visualization and Visual Analysis for Multiphase Flow. *Proceedings of the DIPSI Workshop 2019,* *Università degli studi di Bergamo,* 25–27, 2019.

[4] Alexander Straub, and Thomas Ertl. Visualization Techniques for Droplet Interfaces and Multiphase Flow. *Droplet Interactions and Spray Processes,* *Springer International Publishing,* 121: 203–214, 2020.

