This code computes the proposed Hamiltonian discretisation described in the paper and demonstrates the orthonormality of the discrete Hamiltonian operators relative to the existing method. Loads different meshes (e.g., `centaur3.off`, `Armadillo.off`) using the `read_off` function developed by Gabriel Peyré. The function is available on the MATLAB File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph/files/toolbox_graph/read_off.m. The script stores the mesh vertices and faces (an example is shown below).

## `computeLBOandR.m`

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

## `plot_mesh_with_potential.m`

This script assigns a scalar potential function to the mesh vertices. Different potential functions can be selected, including **step**, **linear**, **strip**, and **random**.

### Inputs

| Parameter | Description |
|-----------|-------------|
| `V` | `n × 3` matrix containing the Cartesian coordinates of the mesh vertices. |
| `F` | `m × 3` triangle connectivity matrix. Each row contains the indices of the three vertices of a triangular face. |
| `pot` | `n × 1` vector containing the scalar potential values assigned to the mesh vertices. |
### Output

| Parameter | Description |
|-----------|-------------|
| |Assigned Potential on the mesh with appropriate colour.|


## `sparsity_plot.m`

This script normalises the computed eigenvectors, constructs the Gram matrix, and generates the sparsity plot to verify the orthogonality of the eigenfunctions.

### Inputs

| Parameter | Description |
|-----------|-------------|
| `eigvec` | Matrix whose columns are the eigenvectors of the Hamiltonian operator computed using MATLAB's `eigs` function. |

### Outputs

| Parameter | Description |
|-----------|-------------|
| `G` | Gram matrix of the normalised eigenvectors. |
| Sparsity plot | Visualisation of the Gram matrix used to verify the orthogonality of the computed eigenfunctions. |

### Method
1. Assembles the global cotangent stiffness matrix W from the "computeLBOandR.m" code
2. Simultaneously assembles the consistent mass matrix and the potential-weighted matrix(R) using a common sparsity pattern, thereby reducing memory allocation and improving computational efficiency.
3. After evaluating all the matrices, one can compute the generalised eigenvalues and eigenvectors of the Hamiltonian operator H_{prop} = (W+R)
4. For orthogonality provided in the figure.png file, one can use the computed eigenvectors $\Psi_{prop}$ and then compute the Grammian matrix $ G_{diag} = \Psi_{prop}^{T} * A * \Psi_{prop}$. 
### Usage
1. Run "assign_potential.m" to define the potential function on the selected mesh.
2. Run "eigenvalue_eigenvector_computation.m" to:
   . construct the Laplace–Beltrami operator,
   . assemble the Hamiltonian operator,
   . compute the eigenvalues and eigenvectors, and   
3. Run "sparsity_plot.m" to verify the orthogonality of the computed eigenfunctions and reproduce Figure - 5,6 & 7 in the paper.
### Repository Structure

| File | Description |
|------|-------------|
| `plot_mesh_with_potential.m` |Loads different meshes (e.g., `centaur3.off`, `Armadillo.off`) using the `read_off` function developed by Gabriel Peyré. The function is available on the MATLAB File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph/files/toolbox_graph/read_off.m. The script stores the mesh vertices and faces (an example is shown below). Assigns the potential function on the selected mesh. |
| `computeLBOandR.m` | Computes the Laplace–Beltrami operator and the \(R\) matrices. |
| `sparsity_plot.m` |Computes the eigenvalues and eigenvectors of the Hamiltonian operator using MATLAB's `eigs` function and stores the computed eigenvalues and eigenvectors (an example is shown below). Verifies the orthogonality of the computed eigenfunctions obtained using different methods.|

### To reproduce Figure 5 presented in the paper, run the following MATLAB commands.

```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the mesh by using read_off (e.g., centaur3.off or Armadillo.off).
[v,f] = read_off('Armadillo.off');
vMat = v';
fMat = f';

% Assign the potential function.
y = vMat(:,2);

% Normalize to [0,1].
ymin = min(y);
ymax = max(y);
y_norm = (y - ymin) / (ymax - ymin);

Pot = zeros(size(y));

% Define strip width.
w = 0.08;   % 8% thickness

% Strip centers (bottom, middle, top).
centers = [0.25, 0.5, 0.75];

for c = centers
    Pot(abs(y_norm - c) <= w) = 10;   % mu = 10
end

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

[Vprop, Dprop] = eigs(H, M, 200, 'smallestabs');
[Vdiag, Ddiag] = eigs(H1, M, 200, 'smallestabs');
[Vlump, Dlump] = eigs(H2, Ml, 200, 'smallestabs');

[eval_diag, ind] = sort(diag(Ddiag));
evec_diag = V2(:, ind);

[eval_prop, ind] = sort(diag(Dprop));
evec_prop = Vprop(:, ind);

[eval_lump, ind] = sort(diag(Dlump));
evec_lump = Vlump(:, ind);

% Compute sparsity.
for k = 1:100
    norm_factor = sqrt(evec_prop(:,k)' * A * evec_prop(:,k));
    evec_prop(:,k) = evec_prop(:,k) / norm_factor;
end

for k = 1:100
    norm_factor = sqrt(evec_diag(:,k)' * A * evec_diag(:,k));
    evec_diag(:,k) = evec_diag(:,k) / norm_factor;
end

Ir = evec_prop' * A * evec_prop;
Ia = evec_diag' * A * evec_diag;

error_thresh = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, ...
                1e-5, 1e-4, 1e-3, 1e-2, 1e-1];

nnz_vals = zeros(size(error_thresh));
nnz_vals1 = zeros(size(error_thresh));

for i = 1:length(error_thresh)
    threshold = error_thresh(i);
    nnz_vals(i) = nnz(abs(Ir) >= threshold);
    nnz_vals1(i) = nnz(abs(Ia) >= threshold);
end

% Plotting.
figure;
plot(log10(error_thresh), nnz_vals, 'k-s', ...
     log10(error_thresh), nnz_vals1, 'k-o', ...
     'LineWidth', 1.5);

legend('Proposed', 'Diag');
grid on;
axis square;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
```

Figure - 5,6 & 7 can be reproduced by the above usage

### Expected Runtime
Approximately 2 minutes on a standard desktop computer with MATLAB R2024a.

### Contact
For questions regarding the implementation, or if required to produce any other figure, please contact the corresponding author.
