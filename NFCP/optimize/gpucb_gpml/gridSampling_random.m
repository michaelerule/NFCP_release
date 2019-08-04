function X = gridSampling_random (n, lbounds, ubounds)

dim = length(lbounds);
ranges = ubounds - lbounds;
rangesMatrix = repmat(ranges, n, 1);
lboundsMatrix = repmat(lbounds, n, 1);
X = rand(n, dim) .* rangesMatrix + lboundsMatrix;
