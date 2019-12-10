function [out, outi] = sortvectors(in)
%sorts the columns of a 2-d matrix (in)
%involves some lossy compression
[~, rankvec] = sort(in);
[~, outi] = sortrows(in');
% [~, outi] = sortrows(rankvec');
out = in(:, outi);