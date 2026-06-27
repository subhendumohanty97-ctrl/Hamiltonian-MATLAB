## computeLBOandR

This function assembles the finite element matrices associated with the Laplace–Beltrami operator on a triangular surface mesh using piecewise linear (P1) finite elements. Given the mesh connectivity, vertex coordinates, and a scalar potential defined at the mesh vertices, it computes the cotangent stiffness matrix, the consistent mass matrix, and the potential-weighted matrix.

### Inputs

| Parameter | Description |
|-----------|-------------|
| `F` | `m × 3` triangle connectivity matrix. Each row contains the indices of the three vertices of a triangular face. |
| `V` | `n × 3` matrix containing the Cartesian coordinates of the mesh vertices. |
| `pot` | `n × 1` vector of scalar potential values defined at the mesh vertices. |

### Outputs

| Output | Description |
|--------|-------------|
| `C` | Sparse cotangent stiffness matrix corresponding to the discrete Laplace–Beltrami operator. |
| `M` | Sparse consistent mass matrix obtained using P1 finite elements. |
| `R` | Sparse potential-weighted matrix representing the discretization of the potential term. |

### Method

The function performs the following steps:

1. Computes the edge vectors and area of each triangular element.
2. Evaluates the cotangent weights associated with each triangle.
3. Assembles the global cotangent stiffness matrix.
4. Simultaneously assembles the consistent mass matrix and the potential-weighted matrix using a common sparsity pattern, thereby reducing memory allocation and improving computational efficiency.
5. Constructs the global matrices as MATLAB sparse matrices.

### Usage

```matlab
[C, M, R] = computeLBOandR(F, V, pot);
```

### Applications

The assembled matrices can be used in:
- Laplace–Beltrami eigenvalue problems.
- Hamiltonian operator on triangulated surfaces.
- Spectral geometry and shape analysis.
