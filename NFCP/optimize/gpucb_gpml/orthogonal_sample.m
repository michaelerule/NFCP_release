% orthogonal_sample.m - creates an orthogonal sample of 1:k*m in N dimensions
%
% USAGE:
%
% [S, m] = orthogonal_sample(N,k,n, [DEBUG])
%
% where the parameters are all scalar:
%
% N = number of dimensions
% k = number of "large" subdivisions (subspaces) per dimension
% n = number of samples to place in each subspace
%
% DEBUG: if this parameter exists, the results (projected into the first two dimensions
% by discarding the other dimensions) will be plotted.
%
% Allowed ranges:
%
% N >= 1, N integer
% k >= 1, k integer
% n >= 1, n integer
%
% If k = 1, the algorithm reduces to Latin hypercube sampling.
% If N = 1, the algorithm simply produces a random permutation of 1:k.
%
% Let m = n*k^(N-1) the number of bins for one variable in one subspace.
% The total number of samples is always exactly k*m.
% Each component of a sample can take on values 1, 2, ..., k*m.
%
% S (output) = samples; (k*m)-by-N matrix where each row is an N-tuple of numbers in 1:k*m.
% m (output) = number of bins per parameter in one subspace.
%              (This is always n*k^(N-1), but is provided as output for convenience.)
%
% See e.g. Wikipedia for a description of both orthogonal and Latin hypercube sampling:
% http://en.wikipedia.org/wiki/Latin_hypercube_sampling 
%
% Briefly, this is a bit like designing the first stage of a Sudoku puzzle: each subspace must 
% have exactly the same number of samples, and no two samples may occur on the same hyperplane.
%
% Examples:
%
% N = 2 dimensions, k = 3 subspaces per axis, n = 1 sample per subspace. m will be n*k^(N-1) = 1 * 3^(2-1) = 3.
% [S,m] = orthogonal_sample(2, 3, 1, 'show')
%
% Compare the previous example with this Latin hypercube that has 9 samples in total:
% (We choose 9 because in the previous example, k*m = 3*3 = 9.)
% [S,m] = orthogonal_sample(2, 1, 9, 'show')
%
% JJ 2010-09-23
%
function [S, m] = orthogonal_sample(N,k,n,DEBUG)
	% Parameter parsing.
	%
	if nargin < 3  ||  N < 1  ||  k < 1  ||  n < 1
		fprintf(1,'Invalid parameters. Please see "help orthogonal_sample" for instructions.\n');
		S = [];
		m = 0;
		return;
	end

	if nargin > 3
		debug_vis = 1;
    else
        debug_vis = 0;
	end

	% Bullet-proof the parameters.
	if N - floor(N) ~= 0
		fprintf(1,'orthogonal_sample.m: WARNING: non-integer N = %0.3g; discarding decimal part\n', N);
		N = floor(N);
	end
	if k - floor(k) ~= 0
		fprintf(1,'orthogonal_sample.m: WARNING: non-integer k = %0.3g; discarding decimal part\n', k);
		k = floor(k);
	end
	if n - floor(n) ~= 0
		fprintf(1,'orthogonal_sample.m: WARNING: non-integer n = %0.3g; discarding decimal part\n', n);
		n = floor(n);
	end

	% Discussion.

	% Proof that the following algorithm implements orthogonal sampling:
	%
	% * Orthogonal sampling has two properties: Latin hypercube sampling globally, and equal density in each subspace.
	% * The independent index vector generation for each parameter guarantees the Latin hypercube property: some numbers will have been used, and removed from the index vectors, when the next subspace along the same hyperplane is reached. Thus, the same indices cannot be used again for any such subspace. This process continues until each index has been used exactly once.
	% * The orthogonal sampling property is enforced by the fact that each subspace gets exactly one sample generated in one run of the loop. The total number of samples is, by design, divisible by the number of these subspaces. Therefore, each subspace will have the same sample density. 

	% Run time and memory cost:
	%
	% * Exactly k*m samples will be generated. This can be seen from the fact that there are k*m bins per parameter, and they all get filled by exactly one sample.
	% * Thus, runtime is in O(k*m) = O( k * n*k^(N-1) ) = O( n*k^N ). (This isn't as bad as it looks. All it's saying is that a linear number of bins gets filled. This is much less than the total number of bins (k*m)^N - which is why orthogonal sampling is needed in the first place. Orthogonal sampling gives us a reduction in sample count by the factor (k*m)^(N-1).)
	% * Required memory for the final result is (k*m)*N reals (plus some overhead), where the N comes from the fact that each N-tuple generated has N elements. Note that the index vectors also use up k*m*N reals in total (k*N vectors, each with m elements). Thus the memory cost is 2*k*m*N reals plus overhead.
	% * Note that using a more complicated implementation that frees the elements of the index vectors as they are used up probably wouldn't help with the memory usage, because many vector implementations never decrease their storage space even if elements are deleted.
	% * In other programming languages, one might work around this by using linked lists instead of vectors, and arranging the memory allocations for the elements in a very special way (i.e. such that the last ones in memory always get deleted first). By using a linked list for the results, too, and allocating them over the deleted elements of the index vectors (since they shrink at exactly the same rate the results grow), one might be able to bring down the memory usage to k*m*N plus overhead.
	% * Finally, note that in practical situations N, k and m are usually small, so the factor of 2 doesn't really matter.

	% Algorithm.

	% Find necessary number of bins per subspace so that equal nonzero density is possible.
	% A brief analysis shows that in order to exactly fill up all k*m bins for one variable,
	% we must have k*m = n*k^N, i.e...
	m = n*k^(N-1);

	% Create index vectors for each subspace for each parameter. (There are k*N of these.)
	I = zeros(N,k,m);
	for i = 1:N
		for j = 1:k
			% We create random permutations of 1:m here so that in the sampling loop
			% we may simply pick the first element from each vector. This is conceptually
			% somewhat clearer than the alternative of picking a random element from an ordered vector.
			%
			% Note that the permutation property ensures that each integer in 1:m only appears once in each vector.
			I(i,j,:) = randperm(m);
		end
	end

	% Start with an empty result set. We will place the generated samples here.
	S = [];

	L = k*m; % number of samples still waiting for placement
	Ns = k^N; % number of subspaces in total (cartesian product of k subspaces per axis in N dimensions)
	while L > 0
		% Loop over all subspaces, placing one sample in each.
		for j = 1:Ns % index subspaces linearly
			% Find, in each dimension, which subspace we are in.
			pj = zeros(N); % multi-index (vector containing an index in each dimension) for this subspace
			% Simple example: (N,k,n) = (2,3,1) --> pj = 1 1, 2 1, 3 1, 1 2, 2 2, 3 2, 1 3, 2 3, 3 3
			pj = 1 + mod( floor( (j-1) ./  k.^(0:N-1) ), k );

			% Construct one sample point.
			% Note: we must use a for loop because ":" doesn't allow to specify two indices as the same.
			% FIXME: better ways to do this?
			%
			s = zeros(1,N);
			for i = 1:N
				s(i) = I(i, pj(i), 1); % grab first element in the corresponding index vector
				I(i,pj(i),:) = cat(3, I(i,pj(i),2:end), -1); % remove the used element, bump others toward the front
			end

			% Now s contains a sample from (1:m)^N. By its construction,
			% the sample conforms to the Latin hypercube requirement.

			% Compute base index along each dimension: the first element of the current subspace
			% is at the multi-index (a+1) (where the +1 is scalar, added to each element of a)
			%
			a = (pj-1)*m;

			% Add the new sample to the result set.
			S = [ S ; a+s ];
		end

		% Note that we placed exactly Ns samples during the for loop.
		L = L - Ns;
	end

	% Result visualization (for debug and illustrative purposes)
	%
	if debug_vis  &&  N > 1
		figure(1);
		clf;

		% Plot the picked points
		plot(S(:,1), S(:,2),'bo');
		axmax = k*m + 1;

		hold on;

		% Mark the subspaces onto the figure
		for j = 1:k-1
			xy = 0.5 + j*m;
			plot( [xy, xy], [0.5, axmax - 0.5], 'k', 'LineWidth', 2.0 );
			plot( [0.5, axmax - 0.5], [xy, xy], 'k', 'LineWidth', 2.0 );
		end

		% Mark bins
		for j = 1:(k*m)-1
			xy = 0.5 + j;
			plot( [xy, xy], [0.5, axmax - 0.5], 'k');
			plot( [0.5, axmax - 0.5], [xy, xy], 'k');
		end

		% Make a box around the area
		plot( [0.5,         axmax - 0.5], [0.5,         0.5],         'k', 'LineWidth', 2.0  );
		plot( [0.5,         axmax - 0.5], [axmax - 0.5, axmax - 0.5], 'k', 'LineWidth', 2.0  );
		plot( [0.5,         0.5],         [0.5,         axmax - 0.5], 'k', 'LineWidth', 2.0  );
		plot( [axmax - 0.5, axmax - 0.5], [0.5,         axmax - 0.5], 'k', 'LineWidth', 2.0  );

		hold off;

		% Set the axes so that the extreme indices just fit into the view
		axis( [0.5 axmax-0.5 0.5 axmax-0.5 ] );

	end
end