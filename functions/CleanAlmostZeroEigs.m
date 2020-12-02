function out = CleanAlmostZeroEigs(matrix,tol)
[eigvecs, eigvals] = eig(matrix);
neweigvals = diag(eigvals);
sum_before = sum(neweigvals(:));
neweigvals(abs(neweigvals)<tol)=0; % clean small eigenvalues
sum_after = sum(neweigvals(:));
diff = sum_before - sum_after;
idx =find(neweigvals==max(neweigvals));
% if I'm cleaning almost 0 eigenvalues, I would like to preserve the 
% trace of the operator so I'm adding the small things I'm removing
neweigvals(idx) = neweigvals(idx) + diff; 
out = 0;
for i=1:size(neweigvals,1)
   out = out + neweigvals(i) * eigvecs(:,i)*eigvecs(:,i)'; 
end
end

