function out = softhresh(arg, threshold)
%Returns the soft thresholding of matrix arg

out = wthresh(arg,'s',threshold);
end