%{

Use Matlab's minimize functionality to optimize objective function

Might also want to optimize these? 

'reg_state_var',1e-5     ,'Regularize orthogonal to pop. size & spatial correlation'
'reg_count_var',1e-8     ,'Regularization orthogonal to spatial correlations'
'reg_diag'     ,1e-8     ,'Diagonal regularization'
'reg_inverse'  ,1e-5     ,'Diagonal regularization for covariance inverse'
'llreg'        ,0.0      ,'No regularization for log-Determinant in likelihood'
'minrate'      ,1e-4     ,'Minimum intensity in measurement update'


%}



model.update     = 'Laplace';
model.likemethod = 'ML'; 
model.likemethod = 'ELBO'
model.likemethod = 'Laplace'

model     = initializeModel(model);
modelNLL  = @(model) -modelLikelihood(model,ini,xydata);
par2mod   = @(Re,Ra,Rr) initializeModel(applyOptions(model,{'rAA',Re,'linearRates',[0 Ra Rr]},false));
%par2mod   = @(Re,Ra) initializeModel(applyOptions(model,{'rAA',Re,'linearRates',[0 Ra model.linearRates(3)]},false));
obj       = @(pars) modelNLL(par2mod(pars{:}))
coord2par = @exp
par2coord = @log
objective = @(pars) obj(num2cell(coord2par(pars)))

%x0 = [model.rAA model.linearRates(2)];% model.linearRates(3)];
x0 = [model.rAA model.linearRates(2) model.linearRates(3)];
x  = fminsearch(objective,par2coord(x0))

x 

y = coord2par(x)

pars = num2cell(coord2par(x))
newmodel = par2mod(pars{:})

[llsum,infstate,margvar,infe] = stateInfer(ini,newmodel,xydata,simulatedM,{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,25
    'showduration' ,1000  
    'showmaxy'     ,density+5   
    'ratescale'    ,50
    'peakactivity' ,false
    'points'       ,false
    });


%{
[llsum,infstate,margvar,infe] = stateInfer(ini,newmodel,xydata,simulatedM,{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,100
    'showduration' ,1000  
    'showmaxy'     ,density+5   
    'ratescale'    ,500
    'peakactivity' ,false
    'points'       ,false
    });

ii = cell2mat(infstate)
%}
