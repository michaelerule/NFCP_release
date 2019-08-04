% Explore likelihood surface of rate parameters
figure(2)
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model.likemethod = 'ML'; % maximum likelihood approach discards the prior
K            = 1; % lenght of likelihood (debug mode can return separate contributions to likelihood)
names        = {'lPyaQ'};
domean       = false;
dosum        = true;
NGRID        = 15;
model.safety = 2;
model        = initializeModel(model);
objective    = @(model) modelLikelihood(model,ini,xydata);
figure(1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinate sweep the excitability parameter QA -> AA
sweep_Re     = 10.^linspace(-4,-0.5,NGRID);
param2model  = @(Re) initializeModel(applyOptions(model,{'rAA',Re}));
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
if dosum,
    sumll2 = sum(ll2,2)';
    semilogx(sweep_Re,sumll2);
end
if domean,
    meanll2 = mean(ll2,2)';
    semilogx(sweep_Re,meanll2);
    sumll2 = meanll2./K;
end
legend(names{:});
x = model.rAA;
semilogx([x x],ylim(),'r');
[maxll2,maxll2idx] = max(sumll2);
bestRe = sweep_Re(maxll2idx)
semilogx([bestRe bestRe],ylim(),'b');
title('Excitation strength')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% What happens if we filter at these parameters?
%{

METHOD = 'ELBO';%'Laplace';
model.likemethod = METHOD;
newmodel.likemethod = METHOD;

newmodel = initializeModel(applyOptions(model,{'rAA',bestRe},false));
newmodel = initializeModel(newmodel);
[llsum,infstate,margvar,infe] = stateInfer(ini,newmodel,xydata,false,{...
    'doplot'       ,false 
    'iniP'         ,iniP
    });
llsum
[llsum,infstate,margvar,infe] = stateInfer(ini,model,xydata,false,{...
    'doplot'       ,false 
    'iniP'         ,iniP
    });
llsum


model.likemethod = 'ML';
newmodel.likemethod = 'ML';

infA = cell2mat(cellfun(@(M) model.getA'*M,infstate,'UniformOUtput',false));
ntimes = size(infA,2);
%infA = reshape(infA,model.n,model.n,ntimes);

g = model.total_gain;
b = model.adjusted_bias;

clear plus
clear times

x = infA;
x = bsxfun(@times, x, g);
x = bsxfun(@times, x, model.dx*model.dt);
x = bsxfun(@plus,  x, b);

meanx = mean(x,1);

clf
cla
hold all
plot(meanx)

meancounts = cellfun(@(x) size(x,1),xydata)/(model.n.^2);
plot(meancounts)

legend('mean activation','mean counts')


newmodel = initializeModel(applyOptions(model,{'rAA',bestRe},false));
newmodel = initializeModel(newmodel);
[llsum,infstate,margvar,infe] = stateInfer(ini,newmodel,xydata,false,{...
    'doplot'       ,true 
    'upscale'      ,6
    'skipinf'      ,10/dt
    'showduration' ,1800/dt
    'showmaxy'     ,1
    'ratescale'    ,A_color_scale
    'peakactivity' ,false
    'points'       ,false
    'iniP',iniP
    });
llsum

%}



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
if dosum,
    sumlla = sum(lla,2)';
    semilogx(sweep_Ra,sumlla);
end
if domean,
    meanlla = mean(lla,2)';
    semilogx(sweep_Ra,meanlla);
    sumlla = meanlla./K;
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
if dosum,
    sumllr = sum(llr,2)';
    semilogx(sweep_Rr,sumllr);
end
if domean,
    meanllr = mean(llr,2)';
    semilogx(sweep_Rr,meanllr);
    sumllr = meanllr./K;
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
if dosum,
    sumlls = sum(lls,2)';
    semilogx(sweep_s,sumlls);
end
if domean,
    meanlls = mean(lls,2)';
    semilogx(sweep_s,meanlls);
    sumlls = meanlls./K;
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
%{
newmodel = initializeModel(applyOptions(model,{'rAA',bestRe,'sigma',bests,'linearRates',[0 bestRa bestRr]},false));
newmodel = initializeModel(newmodel);
[llsum,infstate,margvar,infe] = stateInfer(ini,newmodel,xydata,false,{...
    'doplot'       ,true 
    'upscale'      ,6
    'skipinf'      ,10/dt
    'showduration' ,1800/dt
    'showmaxy'     ,1
    'ratescale'    ,A_color_scale
    'peakactivity' ,false
    'points'       ,false
    'iniP',iniP
    });
llsum


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare to original; sanity check that likelihood is higher? 
[llsum,infstate,margvar,infe] = stateInfer(ini,model,xydata,false,{...
    'doplot'       ,true 
    'upscale'      ,6
    'skipinf'      ,10/dt
    'showduration' ,1800/dt
    'showmaxy'     ,1
    'ratescale'    ,A_color_scale
    'peakactivity' ,false
    'points'       ,false
    'iniP',iniP
    });
llsum
%}



