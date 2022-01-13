function obj = calcDelayTime(obj)
% calculates the delay time by determining the cross correlation lag
% between the FM reference and instrument frequencies and ROCOFs 

% Get the delay data file
prompt = sprintf('Choose the folder containing the Delay Time data');
obj.getResultsFileList(prompt)

[TimeStamps, RefFreq, RefROCOF, InstrFreq, InstrROCOF] = DlyData(obj);
nP = numel(TimeStamps);  % number of samples


% time vector
t = TimeStamps-TimeStamps(1);

% plot the Frequencies and ROCOFS
figure(obj.fig), obj.fig = obj.fig+1;
% plot frequencies
subplot(2,1,1)
plot(t,RefFreq,t,InstrFreq)
xlim([0,t(end)])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
legend('Ref. Freq','Instr. Freq','Location','northeast') 
set(gca,'FontSize',12)
% plot ROCOFs
subplot(2,1,2)
plot(t,RefROCOF,t,InstrROCOF)
xlim([0,t(end)])
xlabel('Time (s)')
ylabel('ROCOF (Hz)')
legend('Ref. ROCOF','Instr. ROCOF','Location','northeast') 
set(gca,'FontSize',12)


% Frequency Delay
% To perform the cross-correlation, the frequency needs to be normalized around 0.
xc = xcorr(InstrFreq-mean(InstrFreq),RefFreq-mean(RefFreq));
[~,idx]=max(xc);
% U1 = fft(RefFreq-mean(RefFreq));
% U2 = fft(InstrFreq-mean(InstrFreq));
% xc = ifft(U2.*conj(U1));
% [~,idx]=max(xc);



% Get the 3 values to interpolate
R = zeros(1,3);
R(2) = xc(idx);
R(1)=xc(idx-1);
R(3)=xc(idx+1);
k = idx-floor(length(xc)/2);
dly = cosInterp(R,k,nP);
obj.FreqDelay = dly*mean(diff(TimeStamps));

% Plot the frequency cross correlation

%xc = xcorr(InstrFreq-mean(InstrFreq),RefFreq-mean(RefFreq));
K = linspace(-floor(length(xc)/2),floor(length(xc)/2),length(xc));
% [~,idx]=max(xc);
% idx = idx-floor(length(xc)/2);
figure(obj.fig), obj.fig = obj.fig+1;
subplot(2,1,1)
plot (K,xc,'k')
xlim([K(1),K(end)])
xline(k,'r')
xlabel('Lag (in samples or reports)')
ylabel('cross correlation')
title(sprintf('Frequency Delay = %0.4f s',obj.FreqDelay))

% 
% tauF = lagsF*mean(diff(TimeStamps));
% plot(tau,xcF)
% xline(tMax,'r')
% xlabel('Time (s)')
% ylabel('normalized cross correlation')
% title(sprintf('Delay = %1.3f seconds',obj.FreqDelay))
% set(gca,'FontSize',12)

% ROCOF delay
xc = xcorr(InstrROCOF,RefROCOF);
% U1 = fft(RefROCOF);
% U2 = fft(InstrROCOF);
% xc = ifft(U2.*conj(U1));
[~,idx]=max(xc);
% Get the 3 values to interpolate
R = zeros(1,3);
R(2) = xc(idx);
R(1)=xc(idx-1);
R(3)=xc(idx+1);
k = idx-floor(length(xc)/2);
dly = cosInterp(R,k,nP);
%dly = cosInterp(R,idx,nP);
obj.ROCOFDelay = dly*mean(diff(TimeStamps));


% Plot the ROCOF cross correlation
subplot(2,1,2)
plot (K,xc,'k')
xlim([K(1),K(end)])
xline(k,'r')
xlabel('Lag (in samples or reports)')
ylabel('cross correlation')
title(sprintf('ROCOF Delay = %0.4f s',obj.ROCOFDelay))


% shift plots
figure(obj.fig), obj.fig = obj.fig+1;
subplot(2,1,1)
plot(t,RefFreq,t-obj.FreqDelay,InstrFreq)
subplot(2,1,2)
plot(t,RefROCOF,t-obj.ROCOFDelay,InstrROCOF)


% EXPERIMENTATION
% cosine fit interpolation

% % get the max value and the value on either side
% XcF = xcF(idx-1:idx+1);
% wF = acos((XcF(1)+XcF(3))/(2*XcF(2)));
% phiF = atan((XcF(1)+XcF(3))/(2*XcF(2)*sin(wF)));
% idxF = idx - phiF/wF;
% % now interpolate between the two time values
% lagLow = lagsF(floor(idxF));
% lagHi = lagsF(ceil(idxF));
% lagDelta = idxF - floor(idxF);
% td = (lagHi+ lagDelta)* mean(diff(TimeStamps));

end

%% ========================================================================
function [TimeStamps, RefFreq, RefROCOF, InstrFreq, InstrROCOF] = DlyData(obj)
% reads the file and gets the timeseries
C = readcell(cell2mat(obj.dataFiles(1)));

% get all the data
hdr = string(C(1,:));
col = find(hdr=='Timestamp');
TimeStamps = cell2mat(C(2:end,col));
col = find(hdr=='Ref Freq');
RefFreq = cell2mat(C(2:end,col));
col = find(hdr=='Inst. Freq');
InstrFreq = cell2mat(C(2:end,col));
col = find(hdr=='Ref ROCOF');
RefROCOF = cell2mat(C(2:end,col));
col = find(hdr=='Inst ROCOF');
InstrROCOF = cell2mat(C(2:end,col));

end

function delay = cosInterp(R,idx,nP)
omega = acos((R(1)+R(3))/(2*R(2)));
theta = atan((R(1)-R(3))/(2*R(2)*sin(omega)));
c = -theta/omega;
lag = mod(idx-1+floor(nP/2),nP)-floor(nP/2);  % integer part
delay = lag+c;
end

