%% GMRES Part
clear, clc;
% Read matrix A
% Data = load('datas/lns_131.mat');
% A = Data.Problem.A;
% b = sum(A')';
% b = ones(A_rows, 1);
% spy(A);
A = [4 -1 0 0; -1 4 -1 0; 0 -1 4 -1; 0 0 -1 3];
b = [1;1;1;1];
% Initialization
tol = 0.0001;
maxiter = 100;

A_rows = size(A, 1);
X = zeros(A_rows, 1); % initial guess of X = [0 ... 0]

V = CalResidual(A, b, X); % residual vector of the initial guess
normV0 = norm(V);
V(:,1) = V/normV0;

residual = [normV0/norm(b)];

c = zeros(A_rows, 1);
s = zeros(A_rows, 1);

e1 = zeros(A_rows, 1);
e1(1) = 1;
beta = normV0*e1;

% Main GMRES iterations
for iter = 1:maxiter
    [V(:,iter+1), H(1:iter+1, iter)] = Arnoldi(A, V, iter);
    
    [H(1:iter+1, iter), c(iter), s(iter)] = GivensRotate(H(1:iter+1,iter), c, s, iter);
    
    % update the residual
    beta(iter + 1) = -s(iter) * beta(iter);
    beta(iter) = c(iter) * beta(iter);
    res = abs(beta(iter + 1)) / norm(b); % relative residual
    
    residual = [residual; res];
    
    fprintf('Res: %f in %i iter\n', res, iter);
    
    if res <= tol
        break;
    end
end
% compute X
y = H(1:iter, 1:iter) \ beta(1:iter);
X = X + V(:, 1:iter) * y;


%% Function Defination - CalResidual
function [Res] = CalResidual(A, B, X)
Res = abs(B - A*X);
end

%% Function Defination - Arnoldi
function [v, H] = Arnoldi(A, V, iter)
v = A*V(:,iter); % Krylov Vector
for k = 1:iter
    dot_prod = dot(V(:,k), v);
    H(k) = dot_prod;
    v = v - dot_prod*V(:,k);
end
normV = norm(v);
H(iter+1) = normV;
if normV ~= 0
    v = v / H(iter+1);
end
end

%% Function Defination - Back Substitution
function x=BackSub(U,b,n)
x=zeros(n,1);
for j=n:-1:1
    if (U(j,j) ~= 0)
        x(j)=b(j)/U(j,j);
    end
    b(1:j-1)=b(1:j-1)-U(1:j-1,j)*x(j);
end
end

%% Function Defination - Generate Givens Rotation Matrix
function [c, s] = GenGivensMat(v1, v2)
    sq = sqrt(v1^2 + v2^2);
    c = v1 / sq;
    s = v2 / sq;
end

%% Function Defination - Givens Rotate
function [H, c_new, s_new] = GivensRotate(H, c, s, iter)
  % apply old rotationer to the ith column
  for i = 1:iter-1
    tmp =  c(i) * H(i) + s(i) * H(i + 1);
    H(i+1) = -s(i) * H(i) + c(i) * H(i + 1);
    H(i) = tmp;
  end

  % compute and apply new rotation
  [c_new, s_new] = GenGivensMat(H(iter), H(iter + 1));

  % eliminate H(i + 1, i)
  H(iter) = c_new * H(iter) + s_new * H(iter + 1);
  H(iter + 1) = 0.0;
end