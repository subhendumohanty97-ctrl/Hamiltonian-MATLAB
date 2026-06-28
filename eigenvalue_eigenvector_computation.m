close all;
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

% Display the eigenvalues
%disp(diag(D));

% Display the eigenvectors (as columns of V)
%disp(V);
patch(fMat, vMat(:,1), vMat(:,2), vMat(:,3), evecr(:,100), 'EdgeColor', 'k');
axis equal;
xlabel('X');
ylabel('Y');
zlabel('Z');
colorbar;
