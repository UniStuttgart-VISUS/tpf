# STL Reader

Read an STL file (.stl) containing a triangle mesh.

## Parameters

The following parameters can be set by the user:

| Parameter         | Description                                                  | Type    | Accepted values | Default value |
| ----------------- | ------------------------------------------------------------ | ------- | --------------- | ------------- |
| Calculate normals | Calculate normals from the triangle surface, although they are usually provided in the file. | Boolean |                 | False         |

## Output

The output is a triangle mesh, for whose vertices the following data fields are provided:

| Data field     | Description                                                  | Data   | Type       | Remark                                                       |
| -------------- | ------------------------------------------------------------ | ------ | ---------- | ------------------------------------------------------------ |
| Surface normal | The normal defined at each vertex.                           | Vector | Node-based |                                                              |
| Valid          | Use the normal calculation as means to check the validity of the normals provided in the file. | Scalar | Node-based | Only available when compiled with sanity checks and detailed output. |
