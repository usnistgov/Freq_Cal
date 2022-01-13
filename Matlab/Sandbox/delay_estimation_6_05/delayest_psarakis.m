function [ delay] = delayest_psarakis(wl,wr)
% [ delay] = delayest_psarakis(wl,wr)
%Estimates delay using Psarakis and Evangelidis's time delay estimate 
%method.
%wl and wr are the two signals (in row vectors).
%This function uses all elements in both vectors with a circular
%correllation

%Reference
%Psarakis, Emmanouil Z., and Georgios D. Evangelidis. "An enhanced 
%correlation-based method for stereo correspondence with subpixel 
%accuracy." Computer Vision, 2005. ICCV 2005. Tenth IEEE International 
%Conference on. Vol. 1. IEEE, 2005.
%
% Copyright Travis Wiens t.wiens@usask.ca 2015


N=numel(wl);

%normalize and remove mean
wlnorm=(wl-mean(wl))/norm(wl-mean(wl));
wrnorm=(wr-mean(wr))/norm(wr-mean(wr));

rho=ifft(fft(wlnorm).*conj(fft(wrnorm)));%circular xcorr
[~, idx]=max(rho);%find peak

%determine if subsample peak is to left or right of idx
if rho(mod(idx-1-1,N)+1)<rho(mod(idx+1-1,N)+1)
    d=mod(idx+1-1,N)+1;
else
    d=idx;
end
%peak is now between d-1 and d

lambda=1;%using the same window so the ratio of norms is 1 (see paper)
r=wrnorm*circshift(wrnorm,[0 1])';%correlation between signal and one sample shift

%find subsample peak
tau=(rho(mod(d-1-1,N)+1)-r*rho(d))/(lambda*(r*rho(mod(d-1-1,N)+1)-rho(d))+r*rho(d)-rho(mod(d-1-1,N)+1));
%tau=(rho(mod(d-1-1,N)+1)-r*rho(d))/((r-1)*(rho(mod(d-1-1,N)+1)+rho(d)));%equation modified for lambda=1


delay=d-1+tau;



