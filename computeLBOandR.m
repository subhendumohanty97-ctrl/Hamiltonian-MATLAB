function [C, M, R] = computeLBOandR(F,V,pot)
%   C : cotan stiffness matrix
%   M : mass matrix
%   R : potential-weighted triple-product matrix

n = size(V,1);

I = F(:,1); 
J = F(:,2); 
K = F(:,3);

Vi = V(I,:);
Vj = V(J,:);
Vk = V(K,:);

%% ------------------------------------------------------------
% Geometry
%% ------------------------------------------------------------
Eij = Vj - Vi;
Eik = Vk - Vi;
Ejk = Vk - Vj;

cross_ik = cross(Eij,Eik,2);
dblA     = vecnorm(cross_ik,2,2);
A        = 0.5*dblA;

cot_i = dot(Eij,Eik,2)./dblA;
cot_j = dot(-Eij,Ejk,2)./vecnorm(cross(-Eij,Ejk,2),2,2);
cot_k = dot(-Eik,-Ejk,2)./vecnorm(cross(-Eik,-Ejk,2),2,2);

%% ------------------------------------------------------------
% Cotangent matrix (unchanged)
%% ------------------------------------------------------------
rowsC = [I;J;I;K;J;K;I;J;K];
colsC = [J;I;K;I;K;J;I;J;K];
valsC = [-0.5*cot_k;
         -0.5*cot_k;
         -0.5*cot_j;
         -0.5*cot_j;
         -0.5*cot_i;
         -0.5*cot_i;
          0.5*(cot_j+cot_k);
          0.5*(cot_i+cot_k);
          0.5*(cot_i+cot_j)];

C = sparse(rowsC,colsC,valsC,n,n);

%% ------------------------------------------------------------
% FUSED Mass + R assembly
%% ------------------------------------------------------------

% Shared sparsity pattern
rows = [I;I;I;J;J;J;K;K;K];
cols = [I;J;K;I;J;K;I;J;K];

% Mass coefficients
m12 = A/12;

valsM = [
    2*m12;
    m12;
    m12;
    m12;
    2*m12;
    m12;
    m12;
    m12;
    2*m12
];

M = sparse(rows,cols,valsM,n,n);

% Potential values
vI = pot(I);
vJ = pot(J);
vK = pot(K);

% Triple-product coefficients
a10 = A/10;
a30 = A/30;
a60 = A/60;

valsR = [
    a10.*vI + a30.*(vJ+vK);
    a30.*vI + a30.*vJ + a60.*vK;
    a30.*vI + a60.*vJ + a30.*vK;

    a30.*vI + a30.*vJ + a60.*vK;
    a10.*vJ + a30.*(vI+vK);
    a60.*vI + a30.*vJ + a30.*vK;

    a30.*vI + a60.*vJ + a30.*vK;
    a60.*vI + a30.*vJ + a30.*vK;
    a10.*vK + a30.*(vI+vJ)
];

R = sparse(rows,cols,valsR,n,n);

end
