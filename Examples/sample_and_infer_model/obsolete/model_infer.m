%
% (re)Start profiler
profile off; profile clear; profile on;    

model     = initializeModel(model);
objective = @(model) -logLikelihood(model,ini,xydata);


% Multivariate simplex algorithm
fprintf(1,'Performing multivariate simplex search\n');
options = optimset(@fminsearch);
options.TolX = 1e-4;
xtrue  = [model.linearRates' model.rAA model.sigma model.thr model.gain];
x0  = xtrue.*0.0+1e-4;
%x0(4) = 0.5;
fun = @(x) objective(initializeModel(model,{
    'linearRates', x(1:model.ntransitions)
    'rAA'        , x(model.ntransitions+1)
    'sigma'      , x(model.ntransitions+2)
    'thr'        , x(model.ntransitions+3)
    'gain'       , x(model.ntransitions+4)}));
xopt = fminsearch(fun,x0,options);
xopt
x0
xtrue


%{
% Grid search
x1 = 1e-6;
x2 = 0.25/model.dt;
K  = 100
range = linspace(x1,x2,K);
fun = @(x) objective(applyOptions(model,{'rAA',x}));
for i=1:K,
    ll(i) = fun(range(i));
    cla;
    plot(ll)
    hold off; fr=getframe(gcf); clear fr; 
end
%}  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explore each coordinate centered on inferred solution

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the search space bounds
%          Re,   Ra,   Rr,  thr, gain, alpha, bias 
LBs    = [  -6    -6    -6    -6    1   -1      0  ];
Ranges = [   6     6     6     6    4    2    1e4  ];
UBs    = LBs + Ranges;
d = numel(LBs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create and objective function

% Construct updated model description based on new parametrs
param2model = @(Re,Ra,Rr,thr,gain,alpha,bias) ...
    initializeModel(...
    applyOptions(model,...
        {...% updated model parameters
            'linearRates',[0 Ra Rr]
            'gain' , gain
            'alpha', alpha
            'bias' , bias
            'thr'  , thr
            'rAA'  , Re},false));

% Map from search space into parameter space
search2param = @(x) {...
    10.^x(1),10.^x(2),10.^x(3),10.^x(4),10.^x(5),10.^x(6),x(7)};

% Workaround matlab's broken syntax (can't say f(g(x){:}))
apply = @(f,args) f(args{:});

% Objective function
obj = @(x) stateInfer(ini,apply(param2model,search2param(x)),xydata,false,{
    'verbose',false
    'cache',false
    });

model.bias = 0;

% Start at true solution
solution = log10([...
    model.rAA,
    model.linearRates(2),
    model.linearRates(3),
    model.thr,
    model.gain,
    model.alpha]);
solution(7) = model.bias;

% Perform a coordinate-wise grid search centered on the 
% inferred solution
NGRID = 10;

N  = d*NGRID;
ll = zeros(N,1);

parfor_progress(N);
parfor ii=0:d*NGRID-1,
    j = mod(ii,NGRID)+1;
    ipar = floor(ii./NGRID)+1;
    a = LBs(ipar);
    b = UBs(ipar);
    scan = linspace(a,b,NGRID);
    newsol = solution(:);
    newsol(ipar) = scan(j);
    ll(ii+1) = obj(newsol);
    parfor_progress;
end
parfor_progress(0);


% Plot coordinate search results
ll = reshape(ll,d,NGRID);
nsub = ceil(sqrt(d));
clf
for ipar=1:d,
    a = LBs(ipar);
    b = UBs(ipar);
    scan = linspace(a,b,NGRID);
    subplot(nsub,nsub,ipar);
    cla
    plot(scan,ll(ipar,:));
    x = solution(ipar);
    hold on
    plot([x x],ylim());
end

%}
