function ll_total = logLikelihood(model,ini,xydata,varargin)
    %{
    
    Debugging copy: return contributions to log-likelihood for 
    all time points to figure out what is wrong!
    
    Log-likelihood function for neural-field cox-process. Wraps the
    `stateInfer` function. 
    
    Peforms sanity checks on parameters and
    returns a large negative value if the parameters are outside the 
    region for which filtering can accurately estimate the likelihood
    
    Parameters
    ----------
    model : `struct`
        `model` struct; should be pre-initialized. See `initializeModel`.
    ini : `vector`
        Filtering initial conditions for the model state-space mean 
        concentration. (Initial covariance is set heuristically using
        model parameters)
    xydata : `cell`
        Cell array of list of (x,y) points in [0,1]^2; the point-process 
        observations
        
    Other Parameters
    ----------------
    errornll : `float`, default -1e60
        Large negative likelihood value to return in the event that
        filtering cannot be performed using the provided model parameters.
    doplot : `bool`, default `false`
        Option forwarded to `stateInfer`;
        Whether to plot states while filtering
    save_figure : `bool`, default `false`
        Option forwarded to `stateInfer`;
        Whether to save plots (if doplot is on) for animation
    cache :  `bool`, default `false`
        Option forwarded to `stateInfer`;
        Whether to store computed states and likelihood, to be retrieved
        later if the function is called with the same model and dataset.
        to avoid cache collisions, store in the model structure some
        identifying information that indicates how xypoints and ini were
        generated.
    upscale : `int`, default 8
        Option forwarded to `stateInfer`;
        Amount of upsampling to perform on shown arrays, if plotting is on.
    skipinf : `int`, default 10   
        Option forwarded to `stateInfer`;
        Skip every `skipinf` ("skip-inferred") frames if plotting is turned
        on.
    showduration : `int`, default 200 
        Option forwarded to `stateInfer`; 
        Total number of time-steps to plot (if plotting)
    showmaxy : `int`, default 2   
        Option forwarded to `stateInfer`;
        Y-axis maximum for plotting inferred states
    showprog : `int`, default false
        Option forwarded to `stateInfer`;
        Whether to print debugging information for every frame
    normll : `bool`, default true
        Whether to normalize log-likelihood by the number of data
        samples.
    timeout : `float`, default -1
        If this is a positive number, monitor seconds elapsed and terminate
        the algorithm early if the runtime exceeds `timeout` seconds. 
    
    Returns
    -------
    ll_total : `float`
        log-likelihood of data given model (up to a constant factor)
    %}
    
    %|===== PARAM ====|= VAL =|
    options = applyOptions({...
        'doplot'       ,false
        'save_figure'  ,false
        'cache'        ,false
        'upscale'      ,8    
        'skipinf'      ,10   
        'showduration' ,200  
        'showmaxy'     ,55   
        'showprog'     ,false
        'errornll'     ,-1e60
        'normll'       ,true
        'timeout'      ,-1
        }, varargin);
 
    % Safety level (model.safety)
    %   0=unchecked
    %   1=patch
    %   2=patch and warn
    %   3=abort
    %
    % Recommend safety level > 0 when optimizing
    % we want errors to be handled here rather than forward up to 
    % the optimization code (which may not print appropriate debugging
    % information)
    % 
    % Using safety = 3 will cause problematic parameter sets to hard-fail,
    % potentially better than wasting CPU cycles trying to filter using
    % numerically infeasible parameters
    %model.safety = 2;
    
    % Logging level (model.verbosity)
    %   0=nothing
    %   1=limited
    %   2=verbose
    %   3=very verbose
    %model.verbosity = 0;
    
    % Return large negative log-likehood in the event of unphysical 
    % parameters. This value is dithered to avoid singular matrices
    % that can occur in the sparse GP code if too many points have 
    % the same likelihood.
    
    %ll_total = options.errornll.*(1+(rand()-0.5)*0.1);
    
    %ll_total={}
    
    if any([model.linearRates(:)' model.thr model.rAA model.dt]<0), 
        printModel(model);
        fprintf(2,'ll=%0.4e\n',ll_total);
        warning('Rates, thresholds, and time-steps should be non-negative');
        %return;
    end;
    if any(model.dt.*[model.linearRates(:)' model.rAA]>0.5), 
        printModel(model);
        fprintf(2,'ll=%0.4e\n',ll_total);
        warning('All rates should be slow relative to time-step model.dt');
        %return;
    end;
    %{
    try
    %}
        [ll_total,infstate,margvar,infe] = stateInfer_debug(ini,model,xydata,false,options);
        %{
        % override ll, use state only (no dynamics)
        % This is a temporary debugging thing
        % does not support under-dispersion correction
        x = cell2mat(infstate);
        A = max(0,real(x(model.nn+1:model.nn*2,1:end)));
        volume = model.dt.*model.dx.*model.gain;
        rate = bsxfun(@times,bsxfun(@plus,A,model.bias),volume);
        for i=1:length(xydata),
            x = pointsToHistogram(xydata{i},model.n,false,false);
            %x = model.blur(x);
            counts(1:model.nn,i) = x(:);
        end
        % Not quite a likelihood
        ll_total = mean(mean(log(max(1e-9,rate)).*counts - rate));
        %}
    %{
    catch e %e is an MException struct
        printModel(model);
        fprintf(2,'ll=%0.4e\n',ll_total);
        fprintf(2,'Identifier:\n%s\n',e.identifier);
        fprintf(2,'Message:\n%s\n',e.message);
        fprintf(2,'%s\n',getReport(e,'extended'));
        fprintf(2,'Continuing, but returning large negative likelihood\n');
    end
    %}
    %fprintf(1,'ll=%0.4e\n',ll_total);
    printModel(model);
end

