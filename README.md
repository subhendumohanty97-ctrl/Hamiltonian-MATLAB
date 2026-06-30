This repository contains code and data for reproducing results from our paper: A finite element method based discretisation of the Hamiltonian operator on triangular meshes, Shape Modeling International 2026. The main contribution in this repository is the code ```ComputeLBOandR.m``` that computes the Cotan, Mass and R matrix as described in the paper for any set of vertices and faces of a closed, manifold, triangular mesh. Triangular meshes stored in the Object file format (.off) can be read using the ```read_off``` function developed by Gabriel Peyré. The function is available on the MATLAB File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph/files/toolbox_graph/read_off.m. It extracts the mesh vertices and faces. 
Other functions in this repository include ```plot_mesh_with_potential.m``` that can be used to plot a triangular mesh with vertices coloured using the input potential function, and ```sparsity_plot.m``` that can be used to verify the (non-) orthogonality of a set of eigenvectors with respect to the inner product defined by the input Mass matrix. Read below for details.

## `computeLBOandR.m`

This function assembles the finite element matrices associated with the Laplace–Beltrami operator on a triangular surface mesh using piecewise linear (P1) finite elements. Given the mesh connectivity, vertex coordinates, and a scalar potential defined at the mesh vertices, it computes the cotangent stiffness matrix, the consistent mass matrix, and the potential-weighted(R) matrix. These matrices can be used to discretise the Hamiltonian operator on closed, manifold, triangular meshes. For details, please refer to the paper https://raw.githubusercontent.com/subhendumohanty97-ctrl/Hamiltonian-MATLAB/main/SMI2026_Preprint_HamiltonianOp.pdf

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
### Usage: [C,M,R] = computeLBOandR(F,V,pot)

## `plot_mesh_with_potential.m`

This function plots the mesh with the user-given potential via colours. It linearly interpolates the colours between vertices.

### Inputs

| Parameter | Description |
|-----------|-------------|
| `V` | `n × 3` matrix containing the Cartesian coordinates of the mesh vertices. |
| `F` | `m × 3` triangle connectivity matrix. Each row contains the indices of the three vertices of a triangular face. |
| `pot` | `n × 1` vector containing the scalar potential values assigned to the mesh vertices. |
### Output

| Parameter | Description |
|-----------|-------------|
| Figure |Assigned Potential on the mesh with appropriate colour.|
### Usage: plot_mesh_with_potential(V,F,pot)

## `sparsity_plot.m`

This script normalises the computed eigenvectors using the given inner product matrix, constructs the Gram matrix, and generates the sparsity plot to verify the orthogonality of the eigenfunctions.

### Inputs

| Parameter | Description |
|-----------|-------------|
| `eigvec` | Matrix whose columns are the eigenvectors of the Hamiltonian operator computed using MATLAB's `eigs` function. |
|`M`|Mass matrix. Given vectors $u$ and $v$, the inner product used is $u^T M v$ |

### Outputs

| Parameter | Description |
|-----------|-------------|
| `G` | Gram matrix of the normalised eigenvectors. |
| Figure (Sparsity plot) | Visualisation of the Gram matrix used to verify the orthogonality of the computed eigenfunctions. |
### Usage: G = sparsity_plot(eigvec,M)

### Method
1. Read a mesh file (all codes assumes closed, triangular, manifold meshes) to obtain vertices and faces in matrices $V$ and $F$, of sizes $n \times 3$ and $m \times 3$, resp.
2. Design a potential function $pot$, a vector of size $n \times 1$. 
3. Compute cotan, mass and R matrix using: `[C,M,R] = computeLBOandR(F,V,pot);`
4. The diagonal entries of the lumped mass matrix can be computed by `M_l = sum(M,2);`
5. The proposed, diag and lumped Hamiltonian discretizations can be computed by `Hprop = W+R; Hdiag=W+M*diag(pot); Hlump = W + diag(M_l.*pot);`
6. You can then compute eigenvalues and eigenvectors of these matrices using the MATLAB function `eig` and work with the corresponding spectral basis. For example use `sparsity_plot` to check orthogonality.
    
### To reproduce Figure 5 presented in the paper, run the following MATLAB commands.

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the mesh by using read_off (e.g., centaur3.off or Armadillo.off).
[v,f] = read_off('Armadillo.off');
vMat = v';
fMat = f';

y = vMat(:,2);

% normalize to [0,1]
ymin = min(y);
ymax = max(y);
y_norm = (y - ymin) / (ymax - ymin);

Pot = zeros(size(y));

% define strip width
w = 0.08;   % 8% thickness

% strip centers (bottom, middle, top)
centers = [0.25, 0.5, 0.75];

for c = centers
    Pot(abs(y_norm - c) <= w) = 10;   % mu = 10
end

plot_mesh_with_potential(vMat, fMat, Pot);

% Compute eigenvalues and eigenvectors.
Hprop = sparse(C + R);
Hdiag = sparse(C + (M .* Pot'));

lumped_mass = sum(M,2);
Hlump = sparse(C + diag(lumped_mass .* Pot));

Ml = diag(lumped_mass);

[Vprop, Dprop] = eigs(Hprop, M, 100, 'smallestabs');
[Vdiag, Ddiag] = eigs(Hdiag, M, 100, 'smallestabs');
[Vlump, Dlump] = eigs(Hlump, Ml, 100, 'smallestabs');

[eval_diag, ind] = sort(diag(Ddiag));
evec_diag = Vdiag(:, ind);

[eval_prop, ind] = sort(diag(Dprop));
evec_prop = Vprop(:, ind);

[eval_lump, ind] = sort(diag(Dlump));
evec_lump = Vlump(:, ind);

% Compute sparsity.
sparsity_plot(evec_prop, evec_diag, M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```
### To reproduce Figure 6 presented in the paper, run the following MATLAB commands.
```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the mesh by using read_off (e.g., centaur3.off or Armadillo.off).
[v,f] = read_off('Armadillo.off');
vMat = v';
fMat = f';

% Assign the potential function.
load('random_and_step_potential.mat');
%This potential is available in the list of files.
Pot = pot;
% Define strip width.
%w = 0.08;   % 8% thickness

% Strip centers (bottom, middle, top).
%centers = [0.25, 0.5, 0.75];

%for c = centers
    %Pot(abs(y_norm - c) <= w) = 10;   % mu = 10
%end

% Compute LBO and R.
[C, M, R] = computeLBOandR(fMat, vMat, Pot);

patch('Faces', fMat, ...
      'Vertices', vMat, ...
      'FaceVertexCData', Pot, ...
      'FaceColor', 'interp', ...
      'EdgeColor', 'k', ...
      'EdgeAlpha', 0.1);

axis square;
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
colorbar;
view(180,270);

% Compute eigenvalues and eigenvectors.
Hprop = sparse(C + R);
Hdiag = sparse(C + (M .* Pot'));

lumped_mass = sum(M,2);
Hlump = sparse(C + diag(lumped_mass .* Pot));

Ml = diag(lumped_mass);

[Vprop, Dprop] = eigs(Hprop, M, 100, 'smallestabs');
[Vdiag, Ddiag] = eigs(Hdiag, M, 100, 'smallestabs');
[Vlump, Dlump] = eigs(Hlump, Ml, 100, 'smallestabs');

[eval_diag, ind] = sort(diag(Ddiag));
evec_diag = Vdiag(:, ind);

[eval_prop, ind] = sort(diag(Dprop));
evec_prop = Vprop(:, ind);

[eval_lump, ind] = sort(diag(Dlump));
evec_lump = Vlump(:, ind);

% Compute sparsity.
sparsity_plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```

### Acknowledgement
The "centaur" mesh data is from the TOSCA dataset (Bronstein, A. M., Bronstein, M. M., & Kimmel, R. (2008). Numerical Geometry of Non-Rigid Shapes. Springer.), "Armadillo" mesh from the Stanford 3D scanning repository (Krishnamurthy, V., & Levoy, M. (1996). Fitting smooth surfaces to dense polygon meshes. In Proceedings of the 23rd annual conference on Computer graphics and interactive techniques (pp. 313-324)).

### Contact
For questions regarding the implementation, or if required to produce any other figure, please contact the corresponding author.
