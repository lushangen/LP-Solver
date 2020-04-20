function[data,info] = LP3033954135(A,b,c)

%Instantiate info and data
info = infoClass;
data = dataClass;

%Setup LP
[m,n] = size(A);
A_hat = [A,eye(m)];
col_A_hat = size(A_hat,2);
B = (n + 1:col_A_hat)';
x_hat = zeros(col_A_hat,1);
x_hat(B) = b;
c_hat = zeros(col_A_hat,1);
c_hat(B) = 1;
T = A_hat(:,B)\[b,eye(m)];
y = T(:,2:end)'*c_hat(B);
T = [T;[c_hat'*x_hat,y']];
flg = 0;

%Check degeneracy
if ~check_degeneracy(A,b)
    disp('This LP is degenerate.');
    info.run = 'Failure';
    info.msg = 'This LP is degenerate.';
    return
end

%Phase I
[T,B,x_hat,~,flg,~,lambda, ITER] = revised_simplex(A_hat,B,c_hat,T,x_hat,y);
info.PhaseI.loop = ITER;
data.PhaseI.obj = c_hat'*x_hat;
data.PhaseI.x = x_hat';

if (flg ~= 0)
    info.run = 'Failure';
    info.msg = 'This LP failed most likely due to degenerate cycling.';
    return;
end

if ~all(ismember(B,1:n))
    disp(['There does not exist a feasible solution due to the following artificial column(s):']);
    disp(B(~ismember(B,1:n)));
    info.cases = 3;
    info.run = 'Failure';
    info.msg = ['There does not exist a feasible solution due to the following artificial column(s):',num2str(B(~ismember(B,1:n)))];
    return;
end

x = x_hat(1:end-m);
T(end,1) = T(1:end-1,1)'*c(B);
y = T(1:end-1,2:end)'*c(B);
T(end,2:end) = y;

%Phase II

[~,B,x,y,flg,t,~,ITER] = revised_simplex(A,B,c,T,x,y);
info.run = 'Success';
info.msg = 'No issues';
info.PhaseII.loop = ITER;
data.PhaseII.Primalobj = c'*x;
data.PhaseII.Dualobj = y'*b;
data.PhaseII.x = x;
data.PhaseII.y = y;
data.PhaseII.z = c-A'*y;

if(flg == 0)
    disp('The optimal primal is x:');
    disp(x);
    disp('The optimal dual is y:');
    disp(y);
    info.cases = 1;
elseif (flg == 2)
    disp('The basic feasible solution is x:');
    disp(x);
    disp('The search direction is t:');
    disp(t);
    disp('The feasible solution x-lambda*t')
    temp = zeros(size(x,1),1);
    temp(B) = t;
    disp(x-lambda*temp);
    info.cases = 2;
    info.msg = 'There is not optimal solution';
    data.PhaseII.t = temp;
end
end

%% Helpers
function [T, B, x, y, flg, t, lam, ITER] = revised_simplex(A, B, c, T, x, y)
simplex = true; 
ITER = 0;
flg = 0;
lam = -1;
t = -1;
while simplex
    y = T(end,2:end)';
    %[zmin,s] = min(c-A'*y);
    obj = y'*A-c';
    s = find(obj > 1e-6);
    s = s(find(~ismember(s,B),1));
    z = obj(s);
    if(isempty(z))
        z = 0;
    end
    
    x = zeros(size(A,2),1);
    x(B)= T(1:end-1,1);
    
    if (abs(z) < 1e-14)
        disp('Simplex Method has converged');
        simplex = false;
        disp(['Optimal Objective = ', num2str(obj),' after ', num2str(ITER), ' iterations']);
       	disp('Displaying Optimal solution x, c-A^T*y and their componentwise product');
        disp([x c-A'*y x.*(c-A'*y)]);
        continue;
    end
    
    t = T(1:end-1,2:end)*A(:,s);
    if (t <= 0)
        disp('This LP has no optimal solution');
        simplex = false;
        flg = 2;
        continue;
    end
    
    [lam, r] = Revisedgetr(size(A,2),B,T,t); 
    
    if (lam == 1)
       disp('This LP is degenerate');
       simplex = 0;
       continue;
    end
    if (r < 1)
        disp('This LP has no lower bound');
        simplex = false;
        continue;
    end
    ITER = ITER + 1;
    f = ['Iteration ', num2str(ITER), ' Obj ', num2str(c'*x), '. Smallest component in c-A^T*y: ', ... 
         num2str(z), ' at s =', num2str(s), '. Component r = ', num2str(r), ' out of basis'];
    disp(f)
    [T,B,flg] = RevisedSimplexTableau(B,r,s,t,z,T);
    if (flg == 1)
       simplex = false;
       continue;
    end
    disp('Current Basis is');
    disp(B');
    disp('Current Tableau is');
    disp(T);
end
end

function [Tnew,Bnew,flg] = RevisedSimplexTableau(B,r,s,t,z,T)
%
% This function updates a RevisedSimplexTableau
% 
%
% On input: 
% B: current basis index
% T: current RevisedTableau
% r: B(r)-th column is to LEAVE the index.
% t: Pivot column
% s: s-th  column is to JOIN the index.
% zmin: s-th component in c-A'*y
%

% On output: 
%
% flg:  flg == 0 indicates SUCCESS in updating,
%       flg == 1 indicates FAILURE in updating,
% Bnew: New basis index
% Tnew: New Tableau
%

%
% initialize flg.
%
flg     = 0;
%
% find dimensions of T.
%
[mt,nt] = size(T); 
%
% Set up Bnew
%
B     = B(:);
Bnew  = [B(1:r-1);s;B(r+1:mt-1)];
%
% Setup Tnew
%
Tnew         = zeros(mt,nt); 
if (t(r) == 0)
%
% This is indication of degeneracy. Quit.
%
    flg = 1;
    return;
end
%
% This is the normal case. Proceed. 
%
Temp              = T(r,:)/t(r);
Tnew(1:mt-1,:)    = T(1:mt-1,:) - t*Temp;
Tnew(r,:)         = Temp;
theta0            = z/t(r);
Tnew(mt,:)        = T(mt,:) - theta0*T(r,:);
end

function [lam, r] = Revisedgetr(n,B,T,t)
%
% find the index to kick out
%
% On input: 
% B: current basis index
% T: current Revised Tableau
% t: current pivot column
% s: s-th  column is to JOIN the index.
% n: number of unknowns.
%

% On output: 
%
% r: B(r)-th column is to LEAVE the index.
%    r < 0 indicates unbounded LP.
%

x   = zeros(n,1);
x(B)= T(1:end-1,1);
if (max(t)<n*eps)
    r = -1;
    return;
end
mask        = find(t>0);
[lam, r] = min(x(B(mask))./t(mask));
r           = mask(r);
end

function result = check_degeneracy(A,b)
result = true;
if (rank(A) ~= size(A,1))
    result = false;
    return;
end
for lc = 1:size(A,1)
    nck = combnk(1:size(A,2),lc);
    for i = 1:size(nck,1)
        if (norm(A(:,nck(i))*(A(:,nck(i))\b) - b) < 1e-14)
            result = false;
            return;
        end
    end
end
end
