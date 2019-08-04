function ll_total = modelLikelihood(model,ini,xydata,varargin)
    %{
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
        'errorll'      ,-1e60
        'normll'       ,true
        'timeout'      ,-1
        'catch_errors' ,true
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
    ll_total = options.errorll.*(1+(rand()-0.5)*0.1);
    if any([model.linearRates(:)' model.thr model.rAA model.dt]<0), 
        printModel(model);
        fprintf(2,'ll=%0.4e\n',ll_total);
        warning('Rates, thresholds, and time-steps should be non-negative');
        return;
    end;
    if any(model.dt.*[model.linearRates(:)' model.rAA]>1.0), 
        printModel(model);
        fprintf(2,'ll=%0.4e\n',ll_total);
        warning('All rates should be slow relative to time-step model.dt');
        return;
    end;
    if options.catch_errors,
        try
            [ll_total,infstate,margvar,infe] = stateInfer(ini,model,xydata,false,options);
        catch e %e is an MException struct
            printModel(model);
            fprintf(2,'ll=%0.4e\n',ll_total);
            fprintf(2,'Identifier:\n%s\n',e.identifier);
            fprintf(2,'Message:\n%s\n',e.message);
            fprintf(2,'%s\n',getReport(e,'extended'));
            fprintf(2,'Continuing, but returning large negative likelihood\n');
        end
    else,
        [ll_total,infstate,margvar,infe] = stateInfer(ini,model,xydata,false,options);
    end
    fprintf(1,'ll=%0.4e\n',ll_total);
    printModel(model);
end

