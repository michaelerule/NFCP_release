% Explore likelihood surface of rate parameters
figure(2)
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heuristic initial conditions: 
% run backwards first to get better initials?
fprintf('Approximating initial conditions...\n');
[llsum,infstate,margvar,infe] = stateInfer(ini,model,fliplr(xydata),false,{'doplot',false});
new_ini = infstate{end-1};
ini = new_ini
vv = margvar{end-1};
iniV = diag(vv);
iniP = cinv(iniV);
model.ini_count_var = 0;
model.ini_state_var = 0;
model.ini_reg_diag  = mean(vv);
model = initializeModel(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HARD CODING FORCING THIS TODO 
model.alpha = 0.8
model = initializeModel(model);

% Gain is drifting too high, allowing saturating of refractory state
% Forbid this!
model.gain = 26400%0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up for parameter inference and plotting
model.likemethod = 'ML'; % maximum likelihood approach discards the prior
K            = 1; % lenght of likelihood (debug mode can return separate contributions to likelihood)
names        = {'lPyaQ'};
NGRID        = 15;
model.safety = 2;
model        = initializeModel(model);
objective    = @(model) modelLikelihood(model,ini,xydata);
figure(1)

% Generic "set the model" function
% params is a cell array of key,value pairs
setparams = @(params) initializeModel(applyOptions(model,params,false));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Something is wrong, the recovery rate should be around 0.01, not its 
%absolute minimum of 1e-6 !
% HARD CODING FOR NOW
model = setparams({'linearRates',[0 model.linearRates(2) 0.01]});

% DISABLED: without a prior, the model will try to make the gain extremely high
% Adjust the gain parameter to better match observations
% [llsum,infstate,margvar,infe] = stateInfer(new_ini,model,xydata,false,{'doplot',false,'iniP',iniP});
% model = repairGain(model,infstate,xydata);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate sweep the excitability parameter QA -> AA
sweep_Re     = 10.^linspace(-4,-0.5,NGRID);
param2model  = @(Re) setparams({'rAA',Re});
ll2          = {}
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
semilogx(sweep_Re,sumll2);
legend(names{:});
x = model.rAA;
semilogx([x x],ylim(),'r');
[maxll2,maxll2idx] = max(sumll2);
bestRe = sweep_Re(maxll2idx)
semilogx([bestRe bestRe],ylim(),'b');
title('Excitation strength')
model = param2model(bestRe);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate sweep the recovery parameter A->R
param2model = @(Ra) setparams({'linearRates',[0 Ra model.linearRates(3)]});
sweep_Ra    = 10.^linspace(-4,0,NGRID);
lla         = {};
parfor i=1:NGRID,
    lla{i} = objective(param2model(sweep_Ra(i)))
end
lla = reshape(cell2mat(lla),K,NGRID)';
figure(2); subplot(222); hold off; cla; 
semilogx(repmat(sweep_Ra,[K 1])',lla); hold on;
sumlla = sum(lla,2)';
semilogx(sweep_Ra,sumlla);
legend(names{:});
x = model.linearRates(2);
semilogx([x x],ylim(),'r');
[maxlla,maxllaidx] = max(sumlla);
bestRa = sweep_Ra(maxllaidx)
semilogx([bestRa bestRa],ylim(),'b');
title('Active -> Refractory')
model = param2model(bestRa);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate sweep interaction radius
param2model = @(s) setparams({'sigma',s});
sweep_s     = 10.^linspace(-2,0,NGRID);
lls         = {};
parfor i=1:NGRID,
    lls{i} = objective(param2model(sweep_s(i)))
end
lls = reshape(cell2mat(lls),K,NGRID)';
figure(2); subplot(224); hold off; cla; 
semilogx(repmat(sweep_s,[K 1])',lls); hold on;
sumlls = sum(lls,2)';
semilogx(sweep_s,sumlls);
legend(names{:});
x = model.linearRates(2);
semilogx([x x],ylim(),'r');
[maxlls,maxllsidx] = max(sumlls);
bests = sweep_s(maxllsidx)
semilogx([bests bests],ylim(),'b');
title('Interaction radius')
model = param2model(bests);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate sweep the R->Q rate
% Disabled for now: REPLACED WITH ALPHA SWEEP
%{
param2model = @(Rr) initializeModel(applyOptions(model,{'linearRates',[0 model.linearRates(2) Rr]},false));
sweep_Rr    = 10.^linspace(-6,0,NGRID);
llr         = {};
parfor i=1:NGRID,
    llr{i} = objective(param2model(sweep_Rr(i)))
end
llr = reshape(cell2mat(llr),K,NGRID)';
figure(2); subplot(223); hold off; cla; 
semilogx(repmat(sweep_Rr,[K 1])',llr); hold on
sumllr = sum(llr,2)';
semilogx(sweep_Rr,sumllr);
legend(names{:});
x = model.linearRates(3);
semilogx([x x],ylim(),'r');
[maxllr,maxllridx] = max(sumllr);
bestRr = sweep_Rr(maxllridx)
semilogx([bestRr bestRr],ylim(),'b');
title('Slow recovery')
model = param2model(bestRr);
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate sweep dispersion parameter (alpha)
% THIS REPLACES THE SLOW TIMESCALE SWEEP
% CANNOT USE THIS WITHOUT A PRIOR!

%{
param2model = @(alpha) setparams({'alpha',alpha});
sweep_a     = 2.^linspace(-1,1,NGRID);
lla         = {};
parfor i=1:NGRID,
    lla{i} = objective(param2model(sweep_a(i)))
end
lla = reshape(cell2mat(lla),K,NGRID)';
figure(2); subplot(223); hold off; cla; 
semilogx(repmat(sweep_a,[K 1])',lla); hold on;
sumlla = sum(lla,2)';
semilogx(sweep_a,sumlla);
legend(names{:});
x=model.alpha; semilogx([x x],ylim(),'r');
[maxlla,maxllaidx] = max(sumlla);
besta = sweep_a(maxllaidx)
semilogx([besta besta],ylim(),'b');
title('Dispersion parameter')
model = param2model(besta);
%}


% --------------------------------------------------------------------------
% Adjust the gain parameter to better match observations
%[llsum,infstate,margvar,infe] = stateInfer(new_ini,model,xydata,false,{'doplot',false,'iniP',iniP});
%model = repairGain(model,infstate,xydata);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What happens if we filter at these parameters?
%{
newmodel = initializeModel(applyOptions(model,{'rAA',bestRe,'sigma',bests,'linearRates',[0 bestRa bestRr]},false));
newmodel = initializeModel(newmodel);
ini = new_ini
[llsum,infstate,margvar,infe] = stateInfer(ini,newmodel,xydata,false,{...
    'doplot'       ,true 
    'upscale'      ,6
    'skipinf'      ,4/dt
    'showduration' ,1800/dt
    'showmaxy'     ,1
    'ratescale'    ,A_color_scale
    'peakactivity' ,false
    'points'       ,false
    'iniP',iniP
    });
llsum

%}
