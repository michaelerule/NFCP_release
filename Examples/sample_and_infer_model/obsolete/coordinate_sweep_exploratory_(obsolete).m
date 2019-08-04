% Explore likelihood surface of rate parameters

clf

NGRID = 20

model.update        = 'Laplace-subspace'
model.maxiter       = 10
model.llreg         = 1e-7
model.minrate       = 1e-3
model.reg_inverse   = 1e-6
model.reg_state_var = 1e-6

model         = initializeModel(model);
objective     = @(model) -logLikelihood(model,ini,xydata);

% Coordinate sweep the excitability parameter QA -> AA
sweep_Re     = 10.^linspace(-4,0,NGRID);
param2model  = @(Re) initializeModel(applyOptions(model,{'rAA', Re}));
ll2          = 1./zeros(NGRID,1);
parfor i=1:NGRID,
    ll2(i) = objective(param2model(sweep_Re(i)))
end
subplot(131); hold off; cla; 
semilogx(sweep_Re,ll2);
hold on
x = model.rAA;
semilogx([x x],ylim(),'r');
[minll2,minll2idx] = min(ll2);
bestRe = sweep_Re(minll2idx)
semilogx([bestRe bestRe],ylim(),'b');


% Sanity check the data-only likelihoods
for i=1:length(xydata),
    x = pointsToHistogram(xydata{i},model.n,false,false);
    %x = model.blur(x);
    counts(1:model.nn,i) = x(:);
end

% Ground truth
rate0 = cell2mat(rates);
llx0  = mean(mean(log(max(1e-9,rate0)).*counts - rate0))

% Ground truth sanity check
x      = cell2mat(simulatedM);
A0b    = max(0,real(x(model.nn+1:model.nn*2,1:end)));
rate0b = bsxfun(@times,bsxfun(@plus,A0b,model.bias),model.dt.*model.dx.*model.gain);
llx0b  = mean(mean(log(max(1e-9,rate0b)).*counts - rate0b))

% Inferred using true parameters
model.update  = 'Laplace-subspace'; % Measurement update method
model.cutoff  = false;
model         = initializeModel(model);
[llsum,infstate,margvar,infe] = stateInfer(ini,model,xydata,simulatedM,{...
    'doplot'       ,false 
    });
x     = cell2mat(infstate);
A1    = max(0,real(x(model.nn+1:model.nn*2,1:end)));
rate1 = bsxfun(@times,bsxfun(@plus,A1,model.bias),model.dt.*model.dx.*model.gain);
llx1  = mean(mean(log(max(1e-9,rate1)).*counts - rate1))

% Inferred using estimated paramters
newmodel         = initializeModel(applyOptions(model,{'rAA',bestRe},false));
newmodel.update  = 'Laplace-subspace'; % Measurement update method
newmodel.cutoff  = false;
newmodel         = initializeModel(newmodel);
[llsum2,infstate2,margvar2,infe2] = stateInfer(ini,newmodel,xydata,simulatedM,{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,100
    'showduration' ,1000  
    'showmaxy'     ,density+5   
    'ratescale'    ,50
    'peakactivity' ,false
    'points'       ,false
    });
x     = cell2mat(infstate2);
A2    = max(0,real(x(newmodel.nn+1:newmodel.nn*2,1:end)));
rate2 = bsxfun(@times,bsxfun(@plus,A2,newmodel.bias),newmodel.dt.*newmodel.dx.*newmodel.gain);
llx2  = mean(mean(log(max(1e-9,rate2)).*counts - rate2))

llx0
llx0b
llx1
llx2

% Coordinate sweep the recovery parameter A->R
param2model = @(Ra) initializeModel(applyOptions(model,{'linearRates',[0 Ra model.linearRates(3)]},false));
sweep_Ra    = 10.^linspace(-4,0,NGRID);
lla         = 1./zeros(NGRID,1);
parfor i=1:NGRID,
    lla(i) = objective(param2model(sweep_Ra(i)))
end
subplot(132); hold off; cla; 
semilogx(sweep_Ra,lla)
hold on
x = model.linearRates(3);
semilogx([x x],ylim(),'r');
[minlla,minllaidx] = min(lla);
bestRa = sweep_Ra(minllaidx)
semilogx([bestRa bestRa],ylim(),'b');

% Coordinate sweep the R -> Q rate
param2model = @(Rr) initializeModel(applyOptions(model,{'linearRates',[0 model.linearRates(2) Rr]},false));
sweep_Rr    = 10.^linspace(-6,0,NGRID);
llr         = 1./zeros(NGRID,1);
parfor i=1:NGRID,
    llr(i) = objective(param2model(sweep_Rr(i)))
end
subplot(133); hold off; cla; 
semilogx(sweep_Rr,llr)
hold on
x = model.linearRates(2);
semilogx([x x],ylim(),'r');
[minllr,minllridx] = min(llr);
bestRr = sweep_Rr(minllridx)
semilogx([bestRr bestRr],ylim(),'b');


% errors are large, what happens if we filter at these parameters?
newmodel = initializeModel(applyOptions(model,{'rAA',bestRe,'linearRates',[0 bestRa bestRr]},false));
newmodel.update         = 'Laplace-subspace'; % Measurement update method
newmodel.cutoff         = false
newmodel = initializeModel(newmodel);
[llsum,infstate,margvar,infe] = stateInfer(ini,newmodel,xydata,simulatedM,{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,100
    'showduration' ,1000  
    'showmaxy'     ,density+5   
    'ratescale'    ,50
    'peakactivity' ,false
    'points'       ,false
    });





% Coordinate sweep interaction radius
NGRID = 20
param2model = @(s) initializeModel(applyOptions(model,{'sigma',s},false));
sweep_s     = 10.^linspace(-2,0,NGRID);
lls         = 1./zeros(NGRID,1);
parfor i=1:NGRID,
    lls(i) = objective(param2model(sweep_s(i)))
end
subplot(133); hold off; cla; 
semilogx(sweep_s,lls)
hold on
x = model.linearRates(2);
semilogx([x x],ylim(),'r');
[minlls,minllsidx] = min(lls);
bests = sweep_s(minllsidx)
semilogx([bests bests],ylim(),'b');




% Adjustments for state-inference
% Turn off the spontaneous exictation (infer as extrinsic noise)
model.rQA            = 0.0;
model.linearRates(1) = 0.0;
model.update         = 'Laplace-subspace'; % Measurement update method
model.cutoff         = false
model = initializeModel(model);
fprintf('Animating...')
tic()
profile off; profile clear; profile on;
[llsum,infstate,margvar,infe] = stateInfer(ini,model,xydata,simulatedM,{...
    'doplot'       ,true 
    'upscale'      ,8
    'skipinf'      ,100
    'showduration' ,1000  
    'showmaxy'     ,density+5   
    'ratescale'    ,50
    'peakactivity' ,false
    'points'       ,false
    });
profile off; profile viewer;
toc()

x = cell2mat(infstate);
A = x(model.nn+1:model.nn*2,1:end);
for i=1:length(xydata),
    x = pointsToHistogram(xydata{i},model.n,false,false);
    x = model.blur(x);
    counts(1:model.nn,i) = x(:);
end

% Get regional observation parameters
b          = model.adjusted_bias(:);
ag         = model.premultiplied_gain(:);
volume     = model.dt.*model.dx;
ab         = b./ag;
ab(ag<1e-9)= 0.0;
rate = bsxfun(@times,bsxfun(@plus,A,ab),ag);

mean(rate(:))
mean(counts(:))



y = cell2mat(rates);
x = cell2mat(infstate);
A = max(0,real(x(model.nn+1:model.nn*2,1:end)));
volume = model.dt.*model.dx.*model.gain;
rate = bsxfun(@times,bsxfun(@plus,A,model.bias),volume);
for i=1:length(xydata),
    x = pointsToHistogram(xydata{i},model.n,false,false);
    x = model.blur(x);
    counts(1:model.nn,i) = x(:);
end
clf
hold on
plot(mean(counts,1))
plot(mean(cell2mat(rates),1))
plot(mean(rate,1))
legend('counts','true','infer')



% Not quite a likelihood

ll = mean(mean(log(max(1e-9,rate)).*counts - rate))


