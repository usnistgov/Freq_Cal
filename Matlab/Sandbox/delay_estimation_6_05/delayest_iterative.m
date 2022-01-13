function [ delay, N_i ] = delayest_iterative(u2,u1,N_i_max,d_tol)
%[ delay, N_i ] = delayest_iterative(u2,u1,N_i_max,d_tol);
%Estimates delay using Brent's method with parabolic interpolation to find
%the maximum of the cross correlation function (interpolated using Nyquist 
%sampling theorem).
%N_i_max is the maximum number of iterations
%d_tol is the tolerance to stop at (samples)
%
% Copyright Travis Wiens t.wiens@usask.ca 2015

if nargin<3
    N_i_max=20;%max iterations
end
if nargin<4
    d_tol=1e-5;%stopping tolerance
end

U1=fft(u1);
U2=fft(u2);
XC=U2.*conj(U1);%circular cross correlation
xc=ifft(XC);
[~, idx]=max(xc);%find peak

%% solve for maximum using Brent's method with parabolic interpolation

%neighbors
ax=idx-1-1;%ax, bx and cx are initial brackets and internal point
bx=idx-1+0;
cx=idx-1+1;

f=@(x) -fourier_series(XC,x);%use negative to find maximum

[delay, ~, N_i]=brentmin(ax,bx,cx,f,d_tol,N_i_max);%solve for minimum of f(x)

end

