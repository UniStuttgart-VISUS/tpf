# Two-phase flow framework

This framework aims at helping the implementation of algorithms in a multiphase flow (originally two-phase flow) setting, where data is provided as volume of fluid (VOF) fields and velocity fields.
It is split into multiple parts, where the main part is the header-only library **TPF**.
Additionally, it provides filters for the visualization framework ParaView by wrapping TPF modules.
This way, all modules can be easily ported to other visualization software by providing a wrapper responsible for directing input and output.

## Header-only library

The library is separated into several categories:

| Namespace                                    | Description                                                  |
| -------------------------------------------- | ------------------------------------------------------------ |
| [algorithm](include/tpf/algorithm/readme.md) | General algorithms that are independent of the types provided by this library. |
| data                                         | Mostly spatial data types, such as grids and trees.          |
| exception                                    | Additional exceptions for internal exception handling.       |
| geometry                                     | Implementation of geometric objects using CGAL.              |
| log                                          | Thread-safe logging implementation.                          |
| math                                         | Basic math for geometry and linear algebra.                  |
| module                                       | Base classes and interface for the implementation of modules. |
| mpi                                          | High-level MPI functions for TPF data types.                 |
| performance                                  | Performance logging capabilities.                            |
| policies                                     | Policies through interface classes.                          |
| stdext                                       | Extension and specializations of classes from the standard template library. |
| traits                                       | Traits for templates.                                        |
| [utility](include/tpf/utility/readme.md)     | Small utility functions for templated operations.            |
| vtk                                          | VTK wrapper for conversion of VTK data types to TPF and vice versa. |

### Modules

The following modules are available, most of which are wrapped in a [ParaView plugin](#paraview-plugin):

| Module                      | Description                                                  |
| --------------------------- | ------------------------------------------------------------ |
| Correct VOF                 | Correct the volume of fluid (VOF) field by removing isolated non-zero cells. |
| Droplets                    | Segmentation of droplets in the VOF field, and calculation of droplet-specific properties. |
| Dynamic Droplets            | Static glyph-based visualization of dynamic droplet information over multiple time steps. |
| Flow Field                  | Stream-, streak- and pathline computation, with emphasis on per-droplet frames of reference. |
| Fluid Position              | Calculate a point per cell as the representative fluid position. |
| FS3D Reader                 | Read temporal FS3D files, containing VOF and other scalar fields, as well as velocity fields. |
| FS3D Writer                 | Allows to write rectilinear grids and their attached fields in the FS3D file format. |
| Interface Bending           | Calculate the bending as the change of interface curvature.  |
| Interface Curvature         | Calculate the curvature at the phase interface.              |
| Interface Deformation Glyph | Provide a glyph that can visualize interface stretching and bending. |
| Interface Gradient          | Calculate the gradient at the phase interface.               |
| Interface Stretching        | Calculate how the interface stretches or contracts.          |
| PLIC                        | Reconstruct the interface between two fluid phases with PLIC. |
| PLIC 3                      | Reconstruct the interface between two fluid phases and a solid phase with PLIC. |
| STL Reader                  | Read an STL file containing a triangle mesh.                 |
| Surface Tension             | Calculate the surface tension force at the phase interface.  |
| Test Data                   | Create simple 3D test geometry, represented as VOF field and velocity field. |

## ParaView plugins

Every module is wrapped as a filter or source for ParaView, with multiple modules bundled into a plugin.

| Plugin                             | Description                                                  |
| ---------------------------------- | ------------------------------------------------------------ |
| [basic](#basic-plugin)             | Basic computations for droplets, e.g., interface gradient and interface curvature, mostly on rectilinear grids. |
| [deformation](#deformation-plugin) | Calculation and visualization of droplet interface deformation from interface stretching and bending. |
| [flow](#flow-plugin)               | Flow visualizations, such as stream-, streak- and pathline computation. |
| [octovis](#octovis-plugin)         | Computation on octrees, for handling data from Octo Tiger, emphasizing on stellar data sets. |

### Basic plugin

This plugin contains several basic data readers, writers, sources, and filters.

#### Data readers

These modules are data readers for loading data from files.

| Reader                                                                            | Description                                                                                   |
|-----------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| [FS3D Reader](src/paraview/basic/modules/fs3d_reader/Readme.md)                   | Read temporal FS3D files, containing VOF and other scalar fields, as well as velocity fields. |
| [STL Reader](src/paraview/basic/modules/stl_reader/Readme.md)                     | Read an STL file containing a triangle mesh.                                                  |

#### Data writers

These modules are data writers for writing data to files.

| Filter                                                                            | Description                                                                           |
|-----------------------------------------------------------------------------------|---------------------------------------------------------------------------------------|
| [FS3D Writer](src/paraview/basic/modules/fs3d_writer/Readme.md)                   | Allows to write rectilinear grids and their attached fields in the FS3D file format.  |

#### Sources

These modules are sources for generating data based on input parameters.

| Source                                                                            | Description                                                                           |
|-----------------------------------------------------------------------------------|---------------------------------------------------------------------------------------|
| [Test Data](src/paraview/basic/modules/test_data/Readme.md)                       | Create simple 3D test geometry, represented as VOF field and velocity field.          |

#### Filters

These modules are filters for manipulating input data, mostly defined on rectilinear grids.

| Filter                                                       | Description                                                  |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| [Correct VOF](src/paraview/basic/modules/correct_vof/Readme.md) | Correct the volume of fluid (VOF) field by removing isolated non-zero cells. |
| [Droplets](src/paraview/basic/modules/droplets/Readme.md)    | Segmentation of droplets in the VOF field, and calculation of droplet-specific properties. |
| [Dynamic Droplets](src/paraview/basic/modules/dynamic_droplets/Readme.md) | Static glyph-based visualization of dynamic droplet information over multiple time steps. |
| [Fluid Position](src/paraview/basic/modules/fluid_position/Readme.md) | Calculate a point per cell as the representative fluid position. |
| [Interface Curvature](src/paraview/basic/modules/interface_curvature/Readme.md) | Calculate the curvature at the phase interface.              |
| [Interface Gradient](src/paraview/basic/modules/interface_gradient/Readme.md) | Calculate the gradient at the phase interface.               |
| [PLIC](src/paraview/basic/modules/plic/Readme.md)            | Reconstruct the interface between two fluid phases with PLIC. |
| [PLIC 3](src/paraview/basic/modules/plic3/Readme.md)         | Reconstruct the interface between two fluid phases and a solid phase with PLIC. |
| [Surface Tension](src/paraview/basic/modules/surface_tension/Readme.md) | Calculate the surface tension force at the phase interface.  |

### Deformation plugin

This plugin contains filters for the calculation and visualization of interface deformation.

#### Filters

| Filter                                                       | Description                                                  |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| [Interface Bending](src/paraview/deformation/modules/interface_bending/Readme.md) | Calculate the change in interface curvature of a droplet.    |
| [Interface Deformation](src/paraview/deformation/modules/interface_deformation/Readme.md) | Calculate the interface deformation as interface bending and stretching. |
| [Interface Deformation Glyph](src/paraview/deformation/modules/interface_deformation_glyph/Readme.md) | Visualize a glyph for interface deformation.                 |
| [Interface Stretching](src/paraview/deformation/modules/interface_stretching/Readme.md) | Calculate the interface stretching of a droplet.             |

### Flow plugin

This plugin contains filters for the computation of droplet-local stream-, streak-, and pathlines.

#### Filters

These modules are filters for computing and visualizing flow and flow properties.

| Filter                                                       | Description                                                  |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| [Flow Field](src/paraview/flow/modules/flow_field/Readme.md) | Stream-, streak- and pathline computation, with emphasis on per-droplet frames of reference. |

### OctoVis plugin

This plugin contains specific filters for binary star mergers.

#### Filters

These modules are filters for processing stellar simulation data from Octo Tiger, represented in octrees.

| Filter                                                       | Description                                                  |
| ------------------------------------------------------------ | ------------------------------------------------------------ |
| [Binary Star](src/paraview/octovis/modules/binary_star/Readme.md) | Classification of binary stars in a stellar merger, defined by physical properties. |
| [Flow Field (Octree)](src/paraview/octovis/modules/flow_field_octree/Readme.md) | Specific version of the [Flow Field](src/paraview/flow/modules/flow_field/Readme.md) (see [Flow plugin](#flow-plugin)) to handle simulation data from a stellar merger. |
| [Seed Octree](src/paraview/octovis/modules/seed_octree/Readme.md) | Use the octree data structure to create a seed.              |

# License

This project is published under the MIT license. In the following you can find a list of contributors in alphabetical order.

## List of contributors

- [Alexander Straub](https://github.com/straubar), University of Stuttgart, Germany  
  (alexander.straub@visus.uni-stuttgart.de)
- [Moritz Heinemann](https://github.com/moritz-h), University of Stuttgart, Germany  
  (moritz.heinemann@visus.uni-stuttgart.de)
- [Grzegorz K. Karch](https://github.com/grzegorz-k-karch), University of Stuttgart, Germany

## MIT license

The MIT License (MIT)

Copyright (c) 2017-2020 University of Stuttgart Visualization Research Center

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
