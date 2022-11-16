function out = svdthresh(arg, threshold)
%Returns the soft thresholding of singular values of matrix arg

[U, S, V] = svd(arg, 'econ');
S_diag = diag(S);
no_S_diag = sum(S_diag > threshold);

if no_S_diag == 0
    out = zeros(size(arg));
else
    out = U(:, 1:no_S_diag) * bsxfun(@times, V(:,1:no_S_diag)', S_diag(1:no_S_diag) - threshold);
end