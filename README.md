## computeLBOandR

This function assembles the finite element matrices associated with the Laplace–Beltrami operator on a triangular surface mesh using piecewise linear (P1) finite elements. Given the mesh connectivity, vertex coordinates, and a scalar potential defined at the mesh vertices, it computes the cotangent stiffness matrix, the consistent mass matrix, and the potential-weighted(R) matrix. These matrices can be used to discretise the Hamiltonian operator on triangular closed-manifold meshes. For details, please refer to the paper https://raw.githubusercontent.com/subhendumohanty97-ctrl/Hamiltonian-MATLAB/main/SMI2026_Preprint_HamiltonianOp.pdf

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
| `R` | Sparse potential-weighted(R) matrix representing the discretisation of the potential term. |

### Method

The function performs the following steps:

1. Assembles the global cotangent stiffness matrix.
2. Simultaneously assembles the consistent mass matrix and the potential-weighted matrix(R) using a common sparsity pattern, thereby reducing memory allocation and improving computational efficiency.
3. After evaluating all the matrices, one can compute the generalised eigenvalues and eigenvectors of the Hamiltonian operator H_{prop} = (W+R)
4. For orthogonality provided in the figure.png file, one can use the computed eigenvectors $\Psi_{prop}$ and then compute the Grammian matrix $\Psi_{prop}^{T} * A * \Psi_{prop}$. 
### Usage

### Repository Structure

| File | Description |
|------|-------------|
| `load_mesh.m` | Loads different meshes (`centaur3.off`, `Armadillo.off`). |
| `assign_potential.m` | Assigns the potential function. |
| `computeLBOandR.m` | Computes the Laplace–Beltrami operator and the \(R\) matrices. |
| `eigenvalue_eigenvector_computation.m` | Computes the eigenvalues and eigenvectors of the Hamiltonian operator. |
| `sparsity_plot.m` | Verifies the orthogonality of eigenfunctions obtained using different methods. |

### To reproduce the numerical results presented in the paper, follow the steps below.

1. Open MATLAB and set the repository folder as the current working directory.
2. Run "loadmesh.m" and select the desired mesh (e.g., centaur3.off or Armadillo.off).
3. Run "assign_potential.m" to define the potential function on the selected mesh.
4. Run "eigenvalue_eigenvector_computation.m" to:
   . construct the Laplace–Beltrami operator,
   . assemble the Hamiltonian operator,
   . compute the eigenvalues and eigenvectors, and   
5. Run "sparsity_plot.m" to verify the orthogonality of the computed eigenfunctions and reproduce Figure - 5,6 & 7 in the paper.

The generated figures correspond to the numerical results presented in the paper.


### Reproducing paper results
Figure - 5,6 & 7 can be reproduced by the above usage

### Expected Runtime
Approximately 2 minutes on a standard desktop computer with MATLAB R2024a.

### Contact
For questions regarding the implementation, or if required to produce any other figure, please contact the corresponding author.
