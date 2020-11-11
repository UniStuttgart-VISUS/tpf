# STL Reader

Read an STL file (.stl) containing a triangle mesh.

## Parameters

| Parameter                 | Description                                                                                       | Type          | Accepted values   | Default value |
|---------------------------|---------------------------------------------------------------------------------------------------|---------------|-------------------|---------------|
| Calculate normals         | Calculate normals from the triangle surface, although they are usually provided in the file.      | Boolean       | True, False       | False         |

## Output

The output is a triangle mesh.

| Data field                | Description                                                                                       | Remark                                                                |
|---------------------------|---------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------|
| Surface normal            | The normal defined at each vertex.                                                                |                                                                       |
| Valid                     | Use the normal calculation as means to check the validity of the normals provided in the file.    | Only available when compiled with sanity checks and detailed output.  |
