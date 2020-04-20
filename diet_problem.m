%% Diet Problem
% By Shangen Lu, SID: 3033954135
% Adapted from "The Cost of Substinence" by George J. Stigler (U. of
% Minnesota), table taken from this paper.
% From the paper "By making linear combinations of various commodities 
% it is possible to construct a composite commodity which is superior 
% in all respects to some survivor", allowing us to eliminate all non-
% double-starred foods from the chart, leaving only 9 entries.

input = readtable('diet_data.xlsx');
A = input{:,3:end}';
b = [3;70;.8;12;5;1.8;2.7;18;75];
c = ones(9,1);
[m,n] = size(A);
matlab_solver = linprog(c,-A,-b,[],[],zeros(n,1));
disp(matlab_solver);

% Add slack variables
c = [c;zeros(n,1)];
A = [A,-eye(n)];
[data,info] = LP3033954135(A,b,c);
input.Commodity(data.PhaseII.x(1:9) > 0);
solution = data.PhaseII.x(1:9);
solution(data.PhaseII.x(1:9) > 0);
disp('The solution is');
disp(solution);
disp('With yearly cost');
disp(c'*data.PhaseII.x*365);

