%{
Explore likelihood surface of rate parameters

This uses grid search to find a local maximum in the log-likelihood
%}

close all
figure(2)
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maximum likelihood approach discards the prior
model.update     = 'Laplace';
model.likemethod = 'ML'; 

% number of terms in likelihood
% (debug modes can return separate contributions to likelihood)
K            = 1; 

% names of terms in likelihood (e.g. prior, likelihood)
names        = {'lPyaQ'};

% if there are multiple terms in the liklelihood, 
% whether to average contributions before displaying
domean       = false;

% if there are multiple terms in the liklelihood, 
% whether to sum contributions before displaying
dosum        = true;

% Resolution at which to sweep parameters
NGRID        = 15;

% Set model to repair and warn if covariance/precision matrices singular
model.safety = 2;

model        = initializeModel(model);
objective    = @(model) modelLikelihood(model,ini,xydata);
figure(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate sweep the excitability parameter QA -> AA
sweep_Re    = 10.^linspace(-4,-0.5,NGRID);
param2model = @(Re) initializeModel(applyOptions(model,{'rAA', Re}));
ll2         = {}
parfor i=1:NGRID,
    ll2{i} = objective(param2model(sweep_Re(i)))
end
ll2raw = ll2;
ll2 = cell2mat(ll2);
ll2 = reshape(ll2,K,NGRID)';
figure(2); subplot(221); hold off; cla; 
semilogx(repmat(sweep_Re,[K 1])',ll2);
hold on
sumll2 = sum(ll2,2)';
if dosum,
    semilogx(sweep_Re,sumll2);
end
if domean,
    semilogx(sweep_Re,sumll2./K);
end
legend(names{:});
x = model.rAA;
semilogx([x x],ylim(),'r');
[maxll2,maxll2idx] = max(sumll2);
bestRe = sweep_Re(maxll2idx)
semilogx([bestRe bestRe],ylim(),'b');
title('Excitation strength')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate sweep the recovery parameter A->R
param2model = @(Ra) initializeModel(applyOptions(model,{'linearRates',[0 Ra model.linearRates(3)]},false));
sweep_Ra    = 10.^linspace(-4,0,NGRID);
lla         = {};
parfor i=1:NGRID,
    lla{i} = objective(param2model(sweep_Ra(i)))
end
lla = reshape(cell2mat(lla),K,NGRID)';
figure(2); subplot(222); hold off; cla; 
semilogx(repmat(sweep_Ra,[K 1])',lla);
hold on
sumlla = sum(lla,2)';
if dosum,
    semilogx(sweep_Ra,sumlla);
end
if domean,
    meanlla = sumlla./K;
    semilogx(sweep_Ra,meanlla);
end
legend(names{:});
x = model.linearRates(2);
semilogx([x x],ylim(),'r');
[maxlla,maxllaidx] = max(sumlla);
bestRa = sweep_Ra(maxllaidx)
semilogx([bestRa bestRa],ylim(),'b');
title('Active -> Refractory')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate sweep the R->Q rate
param2model = @(Rr) initializeModel(applyOptions(model,{'linearRates',[0 model.linearRates(2) Rr]},false));
sweep_Rr    = 10.^linspace(-6,0,NGRID);
llr         = {};
parfor i=1:NGRID,
    llr{i} = objective(param2model(sweep_Rr(i)))
end
llr = reshape(cell2mat(llr),K,NGRID)';
figure(2); subplot(223); hold off; cla; 
semilogx(repmat(sweep_Rr,[K 1])',llr);
hold on
sumllr = sum(llr,2)';
if dosum,
    semilogx(sweep_Rr,sumllr);
end
if domean,
    meanllr = sumllr./K;
    semilogx(sweep_Rr,meanllr);
end
legend(names{:});
x = model.linearRates(3);
semilogx([x x],ylim(),'r');
[maxllr,maxllridx] = max(sumllr);
bestRr = sweep_Rr(maxllridx)
semilogx([bestRr bestRr],ylim(),'b');
title('Slow recovery')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate sweep interaction radius
param2model = @(s) initializeModel(applyOptions(model,{'sigma',s},false));
sweep_s     = 10.^linspace(-2,0,NGRID);
lls         = {};
parfor i=1:NGRID,
    lls{i} = objective(param2model(sweep_s(i)))
end
lls = reshape(cell2mat(lls),K,NGRID)';
figure(2); subplot(224); hold off; cla; 
semilogx(repmat(sweep_s,[K 1])',lls);
hold on
sumlls = sum(lls,2)';
if dosum,
    semilogx(sweep_s,sumlls);
end
if domean,
    meanlls = sumlls./K;
    semilogx(sweep_s,meanlls);
end
legend(names{:});
x = model.linearRates(2);
semilogx([x x],ylim(),'r');
[maxlls,maxllsidx] = max(sumlls);
bests = sweep_s(maxllsidx)
semilogx([bests bests],ylim(),'b');
title('Interaction radius')

figure(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What happens if we filter at these parameters?

newmodel = initializeModel(applyOptions(model,{'rAA',bestRe,'linearRates',[0 bestRa bestRr]},false));
newmodel.update = 'Laplace-subspace'; % Measurement update method
newmodel.cutoff = false
newmodel = initializeModel(newmodel);
subplot(211)
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
fprintf('llsum=%f\n',llsum);

