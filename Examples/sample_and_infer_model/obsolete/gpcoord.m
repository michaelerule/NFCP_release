%{
Experimental bounded Gaussian Process coodinate descent

Gaussian process: 
define a basis
a kernel
a prior
%}

figure(2)
clf

% (re)start profiler
profile off; profile clear; profile on;    

% (re)initialize model
model.verbosity = 0;
model = initializeModel(model);

% minimize negative log-likelihood
objective = @(model) -logLikelihood(model,ini,xydata);

% Set up range
x1 = 1e-10;
x2 = 0.5/model.dt;

% Set up basis
NBASIS = 200;
BASIS  = linspace(x1,x2,NBASIS);

% Define kernel
SIGMA = NBASIS./5;
GKERN = gaussian1DblurOperator(NBASIS,SIGMA,0);

% Measurement error model
MALPHA = 1e-3;

% Prior mean and covariance
M = zeros(NBASIS,1);
C = full(GKERN*GKERN')*1e10+eye(NBASIS);
P = pinv(C);

% Parameter to search over
x1 = 1e-16;
x2 = 0.25/model.dt;
x0 = model.rAA;

% Parameter values for each basis function
range = exp(linspace(log(x1),log(x2),NBASIS));
%range = linspace(x1,x2,NBASIS);

% Objective to minimize for this parameter
fun = @(x) objective(applyOptions(model,{'rAA',x}));

% Iniitial guess
i0 = find(range>=x0)
i = i0(1)

for iter=1:500,
    Pupdate = zeros(NBASIS,NBASIS);
    Mupdate = M;
    % try
    figure(1); clf;
    %y = fun(range(i));
    
    ll = stateInfer(ini,applyOptions(model,{'rAA',range(i)}),xydata,simulatedM,{...
    'doplot'       ,true 
    'save_figure'  ,false
    'upscale'      ,8    
    'skipinf'      ,10
    'showduration' ,200  
    'showmaxy'     ,55   
    'ratescale'    ,5    
    });
    y = -ll;
    hold off; fr=getframe(gcf); clear fr; 
    
        Pupdate(i,i) = 1.0/(MALPHA*(abs(y)+1e-9));
        Mupdate(i) = y;
    %{
    catch
        y = 1e6;
        err = sqrt(diag(pinv(P)))*2.0;
        y = M(i)+err(i);
        Pupdate(i,i) = 1e-8;
        Mupdate(i) = y;
    end
    %}
    PMnew = P*M + Pupdate*Mupdate;
    Pnew  = P + Pupdate;
    Mnew  = Pnew\PMnew;
    P = Pnew;
    M = Mnew;
    figure(2);
    clf;
    hold all;
    err = sqrt(diag(pinv(P)))*2.0;%1.96;
    xx = [1:NBASIS, fliplr(1:NBASIS)];
    fill(xx,[M-err; flipud(M+err)]',...
        1,'facecolor','blue','edgecolor','none','facealpha', 0.3);
    plot(M);
    hold off; fr=getframe(gcf); clear fr; 
    [minval, argmin] = min(M-err);
    %if i==argmin,
    %    break
    %end
    i = argmin;
    range(i);
end


