function [counts] = binCounts(model,xypoints)
    %{
    Bin count data to model grid; 
    
    Note: counts are multiplied by the dispersion correction model.alpha
    
    +--------------------------+-------------------------------------------+
    | Variable                 | Units                                     | 
    +==========================+===========================================+
    | returned: counts         | alpha∙spikes / bin_x²∙bin_t               |
    +--------------------------+-------------------------------------------+
    | model.alpha              | spike mean   /   spike standard dev.      |
    +--------------------------+-------------------------------------------+
    | model.inhomogeneous_gain | unitless                                  |
    +--------------------------+-------------------------------------------+
    | model.adjusted_bias      | unitless                                  |
    +--------------------------+-------------------------------------------+
    | model.adjusted_gain      | alpha∙spikes /   bin_x²∙bin_t             |
    +--------------------------+-------------------------------------------+
    | rate (inside obj. fun.)  | spikes       /   bin_x²∙bin_t             |
    +--------------------------+-------------------------------------------+

    
    Parameters
    ----------
    model : model struct, passes information used for binning spikes. 
        At minimum, the following must be defined:
            model.n: positive integer 
                the simulation grid size
            model.interpolate: boolean; 
                Whether interpolate when binning spikes:
        The following optional arguments are possible:
            model.verbosity: int
                set to a value greater than 2 to print logging info
            model.cutoff: boolean
                whether to blur spikes at interaction scale
            model.alpha: positive number
                whether to rescale spike counts for pseudolikelihood
    xypoints : an NSPIKEX x 2 packed arracy of (x,y) coordinates in [0,1]
    
    Returns
    -------
    counts, binned to the grid size specified by the model
    
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create a histogram of 2D xypoints (the observation counts)
    if numel(xypoints)<=0,
        counts = zeros(model.n.^2,1);
    elseif size(xypoints,2)~=2
        if isfield(model,'verbosity') && model.verbosity>2,
            display('Assuming xypoints provided as histogram)');
            size(xypoints)
        end
        counts = xypoints;
    else
        interp = true;
        if isfield(model,'interpolate'),
            interp = model.interpolate;
        end
        counts = pointsToHistogram(xypoints,model.n,model.interpolate,false);
        if isfield(model,'cutoff') && model.cutoff, 
            counts = model.blur(counts); 
        end
    end
    
    % Rescale counts by alpha (dispersion parameter)
    % (rates should also be multiplied by alpha)
    % this is incorporated into an effective gain parameter
    % model.premultiplied_gain
    % which also multiplies the gain by the spatiotemporal volume
    if isfield(model,'alpha'),
        counts = counts*model.alpha;
    end
    
    counts = counts(:);
