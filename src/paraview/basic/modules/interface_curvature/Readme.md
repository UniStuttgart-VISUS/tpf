# Interface Curvature

Calculate the curvature of the interface using the method by Popinet [1], which uses a height-based approach in the volume of fluid (VOF) field, with a fallback on parabola-fitting in degenerate cases.

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
| Interface positions   | Field containing the interface barycenter of each interface cell. | Vector | Cell-based |        |

## Output

The following data fields are appended to the input grid:

| Data field          | Description                                | Data   | Type       | Remark |
| ------------------- | ------------------------------------------ | ------ | ---------- | ------ |
| Interface curvature | The curvature at the interface barycenter. | Scalar | Cell-based |        |

The grid itself is not modified.

---

[1] Stéphane Popinet. An accurate adaptive solver for surface-tension-driven interfacial flows. *Journal of Computational Physics*, 228(16): 5838–5866, 2009.