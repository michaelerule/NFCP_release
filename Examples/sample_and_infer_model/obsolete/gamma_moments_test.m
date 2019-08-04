%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matching gamma moments to a Gaussian
mean = 0.1
variance = 0.01
x = linspace(0,1,100)
px_gaussian = 1/sqrt(2*pi*variance).*exp(-0.5*(x-mean).^2/variance)
clf
plot(x,px_gaussian)
hold on
beta  = mean/variance
alpha = mean*beta
px_gamma = beta.^alpha/gamma(alpha).*x.^(alpha-1).*exp(-beta*x);
plot(x,px_gamma)
legend('Gaussian','gamma')
psi(alpha)-log(beta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check new function  [f,g,Dh]=objectiveSubspaceGamma(x,M,F,y,reg1,bias,gain,variances)

obj = @objectiveSubspaceGamma;


obj = @objectiveExpectedLinear;

% number of channels and mean rate
N = 10;
p = 0.1;
% intensities (prior)
M = rand(N,1)*p;
% intensities (posterior)
x = rand(N,1)*p;
% gains
g = rand(N,1)*4;
% biases
b = rand(N,1);
% Covariance (made up)
C = diag(M.^2);
P = cinv(C);
F = chol(P);
% Boundary regularization
reg1 = 1e-4; 
v = diag(C);
r = (x + b).*g;    
y = poissrnd(r);
[f0,g0,Dh0]=obj(x,M,F,y,reg1,b,g,v);

% numerically check gradient
delta = 1e-6;
for i=1:N
    x2 = x(:);
    x2(i) = x2(i) + delta;
    [f2,g2,Dh2]=obj(x2,M,F,y,reg1,b,g,v);
    g_numeric(i) = (f2-f0)/delta;
end
g0'
g_numeric

% numerically check hessian diagonal
delta = 1e-6;
for i=1:N
    x2 = x(:);
    x2(i) = x2(i) + delta;
    [f2,g2,Dh2]=obj(x2,M,F,y,reg1,b,g,v);
    h_numeric(i,1:N) = (g2-g0)/delta;
end
Dh0
diag(h_numeric)-diag(P)


%{
% Check log-gamma expectation part of the objective function

% Calculation from objective function
xb = x   + b;
xr = xb .* g;
goodchannels = g>0.0;
sp           = (y>1e-3)&(goodchannels);
tau     = 1./v;
xbv     = xb.*tau;
beta    = xbv./g;
alpha   = xbv.*xb;
meanlnr = psi(alpha)-log(beta);
logsp   = sum(y(sp).*meanlnr(sp))

% Condensed function
flogsp   = @(x) sum(y.*(psi((x+b).^2./v)-log((x+b)./v./g)));
f0 = flogsp(x)

% numerically check gradient
delta = 1e-9;
for i=1:N
    x2 = x(:);
    x2(i) = x2(i) + delta;
    f2 = flogsp(x2);
    g_numeric(i) = (f2-f0)/delta;
end
g_numeric'

% Derivative from objective function
xbv2 = 2*xbv;
psi1 = psi(1,alpha);
xb1  = 1./xb;
dd   = xbv2.*psi1 - xb1;
dd   = y.*dd;
dd(~sp) = 0.0;
dd

% Condensed derivative (check is same)
dlogsp   = @(x) y.*((2.*(x+b)./v).*psi(1,(x+b).^2./v)-1./(x+b));
dlogsp(x)

% check matches numeric gradient
dlogsp(x)
g_numeric'

%}

