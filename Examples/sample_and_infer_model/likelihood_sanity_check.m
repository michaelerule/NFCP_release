%{

Sanity check on likelihood calculation

Inferred activations should predict model log-likelihood


%}

%model.likemethod = 'ML'
%model.likemethod = 'ELL'
%model.likemethod = 'Laplace'
model.likemethod = 'ELBO'

[llsum,infstate,margvar,infe] = stateInfer(ini,model,xydata,simulatedM,{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,200
    'showduration' ,1000  
    'showmaxy'     ,density+5   
    'ratescale'    ,25
    'points'       ,false
    'save_figure'  ,false
    });
    
infstate = cell2mat(infstate);
margvar  = cell2mat(margvar );

llsum1 = llsum

%{
    
    +--------------------------+-------------------------------------------+
    | Variable                 | Units                                     | 
    +==========================+===========================================+
    | model.dt                 | seconds      /   bin_t                    |
    +--------------------------+-------------------------------------------+
    | model.n                  | bin_x        /   array_x                  |
    +--------------------------+-------------------------------------------+
    | model.nn                 | bin_x²       /   array_x²                 |
    +--------------------------+-------------------------------------------+
    | model.dx                 | array_x²     /   bin_x²                   |
    +--------------------------+-------------------------------------------+
    | model.volume             | s∙array_x²   /   bin_x²∙bin_t             |
    +--------------------------+-------------------------------------------+
    | model.bias               | spikes       /   bin_x²∙bin_t             |
    +--------------------------+-------------------------------------------+
    | model.gain               | spikes       /   s∙array_x²               |
    +--------------------------+-------------------------------------------+
    | model.alpha              | spike mean   /   spike standard dev.      |
    +--------------------------+-------------------------------------------+
    | model.sigma              | array_x (1=linear size of array)          |
    +--------------------------+-------------------------------------------+
    | model.linearRates        | 1/neuron/s                                |
    +--------------------------+-------------------------------------------+
    | model.rAA                | 1/neuron²/s (1/fraction²/s if normalized) |
    +--------------------------+-------------------------------------------+
    | model.inhomogeneous_gain | unitless                                  |
    +--------------------------+-------------------------------------------+
    | model.adjusted_bias      | unitless                                  |
    +--------------------------+-------------------------------------------+
    | model.adjusted_gain      | alpha∙spikes /   bin_x²∙bin_t             |
    +--------------------------+-------------------------------------------+
    | rate (inside obj. fun.)  | spikes       /   bin_x²∙bin_t             |
    +--------------------------+-------------------------------------------+

    
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


%}


% Re-compute the log-likelihood of counts from inferred means!
counts = {};
for i=1:size(xydata,2),
    % Create a histogram of 2D xypoints (the observation counts)
    % These are premultiplied by alpha
    counts{i} = binCounts(model,xydata{i});  
end
counts = cell2mat(counts);
size(counts)
% units of counts are
% alpha *
% spikes per bin per timestep
% alpha * spikes / bin_x²∙bin_t

x = model.getA'*infstate;
x = bsxfun(@plus,x,model.adjusted_bias);
x = bsxfun(@times,x,model.adjusted_gain);

% Poisson likelihood
counts = counts(:);
x      = x(:);
ll     = counts.*log(x) - x;


sum(ll)/size(xydata,2)

llsum



