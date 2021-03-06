function [x_max, y_max, A]=crit_interp_p(y,x)
%[x_max y_max A]=crit_interp_p(y,x)
%fits a parabola to three points: y=[y(x(1)) y(x(2)) y(x(3))]. 
%Returns the position (x_max) and value (y_max) of the
%interpolated critical point (peak or trough). 
%If x is omitted, it is assumed to be [-1 0 1]
%y=A(1)*x.^2+A(2)*x+A(3)
%
% Copyright Travis Wiens t.wiens@usask.ca 2015

if nargin<2
    x=[-1 0 1];
end


a = (x(3) * (y(2) - y(1)) + x(2) * (y(1) - y(3)) + x(1) * (y(3) - y(2)));
b = (x(3)*x(3) * (y(1) - y(2)) + x(2)*x(2) * (y(3) - y(1)) + x(1)*x(1) * (y(2) - y(3)));


x_max=-b/(2*a);%critical point position

if nargout>1
    d = (x(1) - x(2)) * (x(1) - x(3)) * (x(2) - x(3));%denominator
    c = (x(2) * x(3) * (x(2) - x(3)) * y(1) + x(3) * x(1) * (x(3) - x(1)) * y(2) + x(1) * x(2) * (x(1) - x(2)) * y(3));
    y_max=c/d-b^2/(4*a*d);%critical point value
end

if nargout>2
    A=[a b c]./d;%parabola coefficients
end
