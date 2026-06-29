This code computes the proposed Hamiltonian discretisation described in the paper and demonstrates the orthonormality of the discrete Hamiltonian operators relative to the existing method. 

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

## assign_potential

In the code, we can assign different potential functions to the mesh vertices, such as step, linear, strip, and random potentials.

### Inputs

| Parameter | Description |
|-----------|-------------|
| `F` | `m × 3` triangle connectivity matrix. Each row contains the indices of the three vertices of a triangular face. |
| `V` | `n × 3` matrix containing the Cartesian coordinates of the mesh vertices. |
| `pot` | `n × 1` vector of scalar potential values defined at the mesh vertices. |


## sparsity_plot 
This will compute the Grammian matrix. Once the eigenvectors are computed by the eigs function, the sparsity plot code will normalise the eigenvectors and will compute the Grammian matrix.

### Method
1. Assembles the global cotangent stiffness matrix W from the "computeLBOandR.m" code
2. Simultaneously assembles the consistent mass matrix and the potential-weighted matrix(R) using a common sparsity pattern, thereby reducing memory allocation and improving computational efficiency.
3. After evaluating all the matrices, one can compute the generalised eigenvalues and eigenvectors of the Hamiltonian operator H_{prop} = (W+R)
4. For orthogonality provided in the figure.png file, one can use the computed eigenvectors $\Psi_{prop}$ and then compute the Grammian matrix $\Psi_{prop}^{T} * A * \Psi_{prop}$. 
### Usage

### Repository Structure

| File | Description |
|------|-------------|
| `loadmesh.m` | Loads different meshes (e.g., `centaur3.off`, `Armadillo.off`) using the `read_off` function developed by Gabriel Peyré. The function is available on the MATLAB File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/5355-toolbox-graph/files/toolbox_graph/read_off.m. The script stores the mesh vertices and faces (an example is shown below). |
| `assign_potential.m` | Assigns the potential function on the selected mesh. |
| `computeLBOandR.m` | Computes the Laplace–Beltrami operator and the \(R\) matrices. |
| `eigenvalue_eigenvector_computation.m` | Computes the eigenvalues and eigenvectors of the Hamiltonian operator using MATLAB's `eigs` function and stores the computed eigenvalues and eigenvectors (an example is shown below). |
| `sparsity_plot.m` | Verifies the orthogonality of the computed eigenfunctions obtained using different methods. |

### To reproduce Figure 5 presented in the paper, one can run the following command in the command window of MATLAB.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Read the mesh by using read_off (e.g., centaur3.off or Armadillo.off).
[v,f] =  read_off('Armadillo.off');
vMat = v';
fMat = f';

%Assign the potential function
y = vMat(:,2);

%normalize to [0,1]
ymin = min(y);
ymax = max(y);
y_norm = (y - ymin) / (ymax - ymin);
Pot = zeros(size(y));
%define strip width
w = 0.08;   % 8% thickness
%strip centers (bottom, middle, top)
centers = [0.25, 0.5, 0.75];

for c = centers
    Pot(abs(y_norm - c) <= w) = 10;   % mu = 10
end

%Compute LBO and R
[W, A, R] = computeLBOandR(fMat,vMat,Pot);
patch('Faces',fMat,'Vertices',vMat,'FaceVertexCData',Pot,'FaceColor','interp','EdgeColor','k',...
    'EdgeAlpha',0.1);
axis square;
xlabel('X');
ylabel('Y');
zlabel('Z');
grid on;
%title('Vertices = 2562');
colorbar;
view(180,270);
%Compute eigenvalues and eigenvectors
H = sparse(W+R);
%H1 = sparse((W+A*diag(Pot)));
%keyboard;
H1 = sparse(W + (A .* Pot'));
%H2 = sparse((W+Al*diag(Pot)));
lumped_mass = sum(A,2);
H2 = sparse((W+diag(lumped_mass.*Pot)));
%k = 3; % Number of eigenvalues/eigenvectors to compute
%options.tol = 1e-6; % Tolerance for convergence
Al = diag(lumped_mass);
[V3, D3] = eigs(H, A, 200, 'smallestabs'); % Compute the k largest eigenvalues and eigenvectors
[V2, D2] = eigs(H1, A, 200, 'smallestabs');
[V4, D4] = eigs(H2, Al, 200, 'smallestabs');
[evala, ind] = sort(diag(D2));
eveca = V2(:, ind);
[evalr, ind] = sort(diag(D3));
evecr = V3(:, ind);
[evall, ind] = sort(diag(D4));
evecl = V4(:, ind);
%Computing sparsity
for k = 1:100
    norm_factor = sqrt(evecr(:,k)' * A * evecr(:,k));
    evecr(:,k) = evecr(:,k) / norm_factor;
end
for k = 1:100
    norm_factor = sqrt(eveca(:,k)' * A * eveca(:,k));
    eveca(:,k) = eveca(:,k) / norm_factor;
end

Ir = evecr' * A * evecr;
Ia = eveca' * A * eveca;

error_thresh = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
nnz_vals = zeros(size(error_thresh));  % To store counts
nnz_vals1 = zeros(size(error_thresh));
for i = 1:length(error_thresh)
    threshold = error_thresh(i);  % DO NOT overwrite variable 'error'
    nnz_vals(i) = nnz(abs(Ir) >= threshold);  % Count how many elements differ by >= threshold
    nnz_vals1(i) = nnz(abs(Ia) >= threshold);
end

% Plotting
figure;
plot(log10(error_thresh), nnz_vals, 'k-s', log10(error_thresh), nnz_vals1, 'k-o','LineWidth',1.5)
legend('Proposed','Diag')
%xlabel('Error threshold')
%ylabel('Number of elements greater than threshold')
%title('Sensitivity of Ir vs Ia to threshold')
grid on;
axis square;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for figure-6 use the following command


3. Run "assign_potential.m" to define the potential function on the selected mesh.
4. Run "eigenvalue_eigenvector_computation.m" to:
   . construct the Laplace–Beltrami operator,
   . assemble the Hamiltonian operator,
   . compute the eigenvalues and eigenvectors, and   
5. Run "sparsity_plot.m" to verify the orthogonality of the computed eigenfunctions and reproduce Figure - 5,6 & 7 in the paper.

The generated figures correspond to the numerical results presented in the paper.


Figure - 5,6 & 7 can be reproduced by the above usage

### Expected Runtime
Approximately 2 minutes on a standard desktop computer with MATLAB R2024a.

### Contact
For questions regarding the implementation, or if required to produce any other figure, please contact the corresponding author.
