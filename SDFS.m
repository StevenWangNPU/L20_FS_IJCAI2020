function [W, obj] = SDFS(X, Y, m, k)
% X: d by n
% Y: n by 1
% Written by Steven Wang, email: zhengwangml@gmail.com
d = size(X,1);
W0 = orth(rand(d,m));
[Sb, Sw, ~, ~] = calculate_L(X',Y);
lambda = trace(W0'*Sb*W0)/trace(W0'*Sw*W0);
obj = [];
for iter = 1:50
obj = [obj;lambda];
% The abs of smallest value ensures step 18
eigvalue = eigs(Sb - lambda*Sw,1,'SA');
if eigvalue>=0
    A = Sb - lambda*Sw;
else
    eta = - eigvalue;
    A = Sb - lambda*Sw + eta*eye(d);
end
W = IPU(A,m,k,d,W0);
lambda = trace(W'*Sb*W)/trace(W'*Sw*W);
W0 = W;
end