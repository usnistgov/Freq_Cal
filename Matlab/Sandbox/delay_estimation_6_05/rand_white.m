function x=rand_white(N)
%x=rand_white(N);
%Generates exaclty white noise with variance 1.
%N must be even.
%
% Copyright Travis Wiens t.wiens@usask.ca 2015

if mod(N,2)==1
    error('N must be even')
end

X=ones(1,N);
phase=2*pi*rand(1,N/2-1);
X(2:N/2)=X(2:N/2).*exp(2i*pi*phase);
X((N/2+2):end)=conj(fliplr(X(2:N/2)));
x=sqrt(N)*real(ifft(X));
