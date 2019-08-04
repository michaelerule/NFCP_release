function [ logDensity ] = gmmLogPDF( x, weights, means, Sigmas )
% gmmLogPDF Gives the log of the probability density of x for a given GMM
% [ logDensity ] = gmmLogPDF( x, weights, means, Sigmas )
%
%	x:       the points at witch the density id evaluated
%	weights: the component weights of the GMM
%	means:   the component means of the GMM
%	Sigmas:  the component (diagonal) covariance matrices of the GMM
%
%	logDensity: the log of the probability density
%
% Dimitrios Milios (dmilios@inf.ed.ac.uk)

	dim = size(x, 2);
	k = length(weights);
	
	if size(means, 1) ~= k
		error(['The mean vectors do not correspond to ' num2str(k) ' components']);
	end
	if size(Sigmas, 3) ~= k
		error(['The covariance matrices do not correspond to ' num2str(k) ' components']);
	end
	if size(means, 2) ~= dim
		error('Incorrect dimension for the mean vectors');
	end
	if size(Sigmas, 1) ~= dim || size(Sigmas, 2) ~= dim
		error('Incorrect dimension for the covariance matrices');
	end
	
	
	density = 0;
	for i = 1:k
		mu = means(i, :);
		C = Sigmas(:, :, i);
		w = weights(i);
		density = density + w * mvnpdf(x, mu, C);
	end
	logDensity = log(density);
end

