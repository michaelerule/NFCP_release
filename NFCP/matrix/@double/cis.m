%z=cis(p, r)
%-----------
%
%Fast sine and cosine for Intel Pentium II+.
%
%  Leutenegger Marcel © 20.6.2004
%  Linux port by Stefan Roth
%
%Input:
% p     Phase [rad]
% r     Magnitude {1}
%
%Output:
% z     r.*complex(cos(p),sin(p))
%
%Alternative:
% [c,s] To get the sine and cosine separately.
%
function [c,s]=cis(p, r)
switch nargin
case 0
   fprintf('\nSine and cosine.\n\n\tLeutenegger Marcel © 20.6.2004\n');
case 1
   c=exp(complex(0,p));
case 2
   c=r.*exp(complex(0,p));
otherwise
   error('Incorrect number of arguments.');
end
if nargout == 2
   s=imag(c);
   c=real(c);
end
