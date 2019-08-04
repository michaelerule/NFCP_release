function [model,scan,ll,best]=coordinateSearch(model,xydata,ini,LB,UB,NGRID,set,varargin)
    %{
    Search log-likelihood over a fixed range for a single model parameter. 
    
    An initial model, observations, and initial conditions should be 
    provided. 
    The lower bound (LB), upper bound (UB) and number of search points
    (NGRID) should be specified. 
    A function `set` that accepts a parameter and returns a key,value pair
    to set in the model should be provided. 
    For example set=@(s){'sigma',s}.
    
    Parameters
    ----------
    model : struct
        pre-initialized model struct
    xydata : cell
        Point data
    ini : matrix
        initial conditions for mean concentratios
    LB : float
        search lower bound
    UB : float
        search upper bound
    NGRID : int
        number of search points 
    set : function
        setter function, e.g. set=@(s){'sigma',s}
        
    Other Parameters
    ----------------
    cache :  `bool`, default `false`
        Option forwarded to `stateInfer` via `modelLikelihood`;
        Whether to store computed states and likelihood, to be retrieved
        later if the function is called with the same model and dataset.
        to avoid cache collisions, store in the model structure some
        identifying information that indicates how xypoints and ini were
        generated.
    timeout : `float`, default 600
        Option forwarded to `modelLikelihood`.
        If this is a positive number, monitor seconds elapsed and terminate
        the algorithm early if the runtime exceeds `timeout` seconds. 
    
    Returns
    -------
    model : struct
        New model set to the best parameters found during search
    scan : matrix
        List of log10 parameter values searched
    ll : matrix
        List of normalized log-likelihood values at each parameter setting
    best : float
        Best parameter value found during search
    %}
    %|===== PARAM ====|= VAL =|
    options = applyOptions({...
        'timeout',600
        'cache',false
        }, varargin);
    
    functional_macros; 
    search2par = @(x) num2cell(10.^x);
    par2search = @(x) any2mat(log10(x));
    par2mod    = @(x) initializeModel(applyOptions(model,set(x)));
    obj        = @(x) modelLikelihood(ini,apply(par2mod,search2par(x)),xydata,options);
    LB         = par2search(LB);
    UB         = par2search(UB);
    scan       = linspace(LB,UB,NGRID);
    ll         = zeros(NGRID,1);
    parfor j=1:NGRID,
        fprintf(1,'Processing job %d of %d\n',j,NGRID);
        ll(j) = obj(scan(j));
    end
    [~,i]  = max(ll);
    best   = cellitem(search2par(scan(i)),1);
    model  = par2mod(best);
    fprintf(1,'Best found is %0.4f\n',best);
    
