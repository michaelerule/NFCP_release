%{

Sanity check on likelihood calculation

Inferred activations should predict model log-likelihood


%}


backmodel = initializeModel(model);
% Reverse the linear model
backmodel.dt = -model.dt
backmodel.rQA            = 0.0;
backmodel.rAA            = 0.0;
backmodel.linearRates(1) = 0.0;
backmodel.thr            = 0.0;

backmodel.ini_reg_diag = 1e1;
backmodel.reg_diag = 1e1;
backmodel.reg_state_var = 1e1;

% Heuristic backward inference to get initial conditions
[llsum,infstate,margvar,infe] = stateInfer(ini,backmodel,flip(xydata),flip(simulatedM),{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,10
    'showduration' ,1000  
    'showmaxy'     ,density+5   
    'ratescale'    ,25
    'points'       ,false
    'save_figure'  ,false
    });

%{

backmodel = initializeModel(model);
% Reverse the linear model
backmodel.description = [ ...
%    Q  A R1  rate
     1 -1  0  2e-1  % spontaneous excitation 
     0  1 -1  4e-1  % slow refractory loop
    -1  0  1  32e-4 % slow refractory recovery
     ];
% Remove the excitation model
backmodel.rQA            = 0.0;
%backmodel.rAA            = 0.0;
backmodel.linearRates(1) = 0.0;
backmodel.thr            = 0.0;
backmodel.ini_reg_diag = 1e4;

% Heuristic backward inference to get initial conditions
[llsum,infstate,margvar,infe] = stateInfer(ini,backmodel,flip(xydata),flip(simulatedM),{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,10
    'showduration' ,1000  
    'showmaxy'     ,density+5   
    'ratescale'    ,25
    'points'       ,false
    'save_figure'  ,false
    });

%}


infstate = cell2mat(infstate);
margvar  = cell2mat(margvar );
model.ini_reg_diag =  mean(margvar(1:end,end))
ini2 = infstate(1:end,end);

[llsum,infstate,margvar,infe] = stateInfer(ini2,model,xydata,simulatedM,{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,10
    'showduration' ,1000  
    'showmaxy'     ,density+5   
    'ratescale'    ,25
    'points'       ,false
    'save_figure'  ,false
    });
