%vMat1 = [10*vMat(:,1),10*vMat(:,2),10*vMat(:,3)];
%vMat = Vn;
%fMat = f;
%[W,A] = LBO(vMat, fMat);
%lumped_mass = sum(A,2);
%Al = diag(lumped_mass);
%n = size(vMat);
%mu = 10;
%Pot = zeros(size(vMat,1),1);
%Pot = 10*rand(size(vMat,1),1);
% Pot = (Pot/100)*80 + 20;
%Pot(vMat(:,2)>=20) = mu;
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
%Pot = k1;
[W, A, R] = computeLBOandR(fMat,vMat,Pot);
%lumped_mass = sum(A,2);
%sigma = 1e-6;   % small shift
%[eigvecs, eigvals] = eigs(A, B, k, sigma);
%H2 = sparse((W+diag(lumped_mass.*Pot)));
%Pot = rand(size(vMat,1),1);
%Pot = zeros(size(vMat,1),1);
%Pot = 100*rand(size(vMat,1),1);
% Pot(vMat(:,1)>=-1) = mu;
% Pot(vMat(:,2)>=-1) = mu;
%Pot(vMat(:,3)>=0.4) = mu;

%vmin = min(vMat(:,2));
%vmax = max(vMat(:,2));
%Pot = (vMat(:,2)+vmin)/(vmax+vmin);
% Pot( ...
%     (vMat(:,2) >= -7 & vMat(:,2) <= -4) | ...
%     (vMat(:,2) >= -1 & vMat(:,2) <= 4)  | ...
%     (vMat(:,2) > 7  & vMat(:,2) < 10) ...
% ) = 10;
%ri = randperm(size(vMat,1));
%pot2 = Pot(ri);
%Pot = pot2;

%[R] = computeR2(fMat, vMat, Pot);
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
