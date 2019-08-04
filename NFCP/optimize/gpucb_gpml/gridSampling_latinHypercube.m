function X = gridSampling_latinHypercube (lbounds, ubounds, subdivisions, samplesPerSubspace)
% The total number os samples will be:
% subdivisions^d * samplesPerSubspace, where d is the dimension

dim = length(lbounds);
n = subdivisions^dim * samplesPerSubspace;

X = orthogonal_sample(dim, subdivisions, samplesPerSubspace);
X = X / n; %normalise to [0,1]

ranges = ubounds - lbounds;
rangesMatrix = repmat(ranges, n, 1);
lboundsMatrix = repmat(lbounds, n, 1);
X = X .* rangesMatrix + lboundsMatrix;
