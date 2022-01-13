function delay=delayest_fft(y,u)
%delay=delayest_fft(y,u);
%Estimates the delay between row vectors y and u.
%Assumes periodic signals and an even fft size
%
% Copyright Travis Wiens t.wiens@usask.ca 2015

Y=fft(y);
U=fft(u);

XC=Y.*conj(U);%periodic cross correlation

xc=ifft(XC);
[~, idx]=max(xc);%find peak of xcorr
d_int=idx-1;%delay (integer)
N_p=numel(y);%number of points
N_2=N_p/2-1;%number of valid points

phi_shift_int=-d_int*pi*2/N_p*(1:(N_2));

W=abs(XC(2:(N_p/2)));%weight

%phi=unwrap(angle(XC(2:(N_p/2))))-phi_shift_int;
phi=mod(angle(XC(2:(N_p/2)))-phi_shift_int+pi,2*pi)-pi;%alternate method (faster)

d_frac=-(sum((1:(N_2)).*phi.*W.^2)/sum(((1:(N_2)).*W).^2))*(N_p/2)/pi;%fractional delay

delay=d_int+d_frac;
