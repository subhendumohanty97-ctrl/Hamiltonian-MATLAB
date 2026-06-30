function plot_mesh_with_potential(vMat, fMat, Pot)

%PLOT_MESH_WITH_POTENTIAL Visualize a triangular mesh with a given potential.
%
% Inputs:
%   vMat - n x 3 matrix of mesh vertices.
%   fMat - m x 3 matrix of triangle connectivity.
%   Pot  - n x 1 vector of potential values at the mesh vertices.
%
% Output:
%   Displays the mesh colored according to the potential values.

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
