function seq=log10seq(n,baseseq)

if ~exist('baseseq','var')
    baseseq = [1 2 3 4 5];
end

s = 1;
seq = [];
while s < n
    seq = [seq s*baseseq];
    s = s*10;
end

seq = seq(seq <= n);