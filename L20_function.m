function W = L20_function(A,k,W)
% Written by Steven Wang
for iter = 1:30
        %update P
        pinvAW = pinv(W'*A*W);
        P = A*W*pinvAW*W'*A;
        diagP = diag(P);
        %update W
        [~,index] = sort(diagP,'descend');
        indexW = index(1:k);
        indexO = index(k+1:end);
        M = A*W;
        MP = M(indexW, :);
        OMP = orth(MP);
        W([indexW],:) = OMP;
        W([indexO],:) =0;   
end
