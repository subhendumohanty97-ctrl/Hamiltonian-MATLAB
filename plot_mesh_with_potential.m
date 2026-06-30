function [Pot, W, A, R] = plot_mesh_with_potential(vMat, fMat)

% Assign a strip potential function on the mesh.

% Extract the y-coordinates of the vertices.
y = vMat(:,2);

% Normalize the y-coordinates to the interval [0,1].
ymin = min(y);
ymax = max(y);
y_norm = (y - ymin) / (ymax - ymin);

% Initialize the potential.
Pot = zeros(size(y));

% Define strip width.
w = 0.08;      % 8% thickness

% Define strip centers (bottom, middle, top).
centers = [0.25, 0.5, 0.75];

% Assign the strip potential.
for c = centers
    Pot(abs(y_norm - c) <= w) = 10;
end

% Compute the Laplace–Beltrami operator and R matrix.
[W, A, R] = computeLBOandR(fMat, vMat, Pot);

% Visualize the potential on the mesh.
figure;
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

end
