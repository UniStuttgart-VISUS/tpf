# Surface Tension

Compute the surface tension force at the fluid interface.

*Note that this is an experimental feature and is probably far from physically correct.*

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
| Interface gradients   | Gradients in interface cells.                                | Vector | Cell-based |        |
| Interface curvature   | Field storing the curvature at the interface barycenter.     | Scalar | Cell-based |        |

## Parameters

The following parameters can be set by the user:

| Parameter   | Description                                                  | Type  | Accepted values | Default value |
| ----------- | ------------------------------------------------------------ | ----- | --------------- | ------------- |
| Coefficient | Surface tension coefficient based on the involved species.   | Float | \>0             | 72.75         |
| Density     | Density of the liquid phase.                                 | Float | \>0             | 0.9982        |
| Time step   | Alternative time step for scaling the surface tension force vectors. If set to zero, use time step derived from the dataset. | Float | &#8805;0        | 0.0           |

## Output

The following data fields are appended to the input grid:

| Data field            | Description                                           | Data   | Type       | Remark |
| --------------------- | ----------------------------------------------------- | ------ | ---------- | ------ |
| Surface tension force | Approximate surface tension force in interface cells. | Vector | Cell-based |        |

The grid itself is not modified.

---

[1] Stéphane Popinet. An accurate adaptive solver for surface-tension-driven interfacial flows. *Journal of Computational Physics*, 228(16): 5838–5866, 2009.