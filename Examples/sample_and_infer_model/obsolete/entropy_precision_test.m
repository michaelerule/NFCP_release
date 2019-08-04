%
%{
The entropy after 1 poisson update on an uninformative improper Gamma 
prior should be

y ln gamma y + (1-y) psi y

The entropy of a univariate gaussian in terms of sigma is

0.5 ln | 2 pi e sigma^2 |

Equating and solving for sigma gives

sigma = sqrt(2 pi /e) gamma(y) exp(y + (1-y) psi(y))

Can this be used to predict the precision update for multivariate gaussian
poisson observations? 
%}
close all; clear all;

y = linspace(0.1,5,100);

h = y + log(gamma(y)) + (1-y).*psi(y);

plot(y,h);
hold on

e = exp(1);

s  =sqrt(2*pi/e).*gamma(y).*exp(y+(1-y).*psi(y));

plot(y,s.^2);

%plot(y,s.^-2)

legend('Gamma entropy','Variance')%,'Precision')o
