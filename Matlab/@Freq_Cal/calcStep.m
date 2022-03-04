function calcStep(obj)
% calculates the settling time using instrument data taken during a frequency step
% The data is from a series of tests with the step time decremented by 1/10 reporting cycle
% This function reads the instrumnt data, interleaves the data and calculated the 10%
% and 90% boundaries and the settling time.

% Get the measuring range data files
prompt = sprintf('Choose the folder containing the Phase Step data');
obj.getResultsFileList(prompt)

% pre-allocate the array using the data from the first file (assume all files are the same number of rows)
C = readcell(cell2mat(obj.dataFiles(1)));
hdr = string(C(1,:));
%col = find(hdr=='FE');
tCol = find(hdr=='Timestamp',1,'first');

% read all the data into a cell array.  some data is longer than others
data = cell(1,numel(obj.dataFiles));
C = readcell(cell2mat(obj.dataFiles(1)));
%data{1,1} = C(2:end,col);
tData{1,1} = C(2:end,tCol);
nData = length(tData{1,1});

for i = 2:numel(obj.dataFiles)
    C = readcell(cell2mat(obj.dataFiles(i)));
    tData{1,i} = C(2:end,tCol); 
    if length(tData{1,i})<nData,nData = length(tData{1,i});end % if the data files are not the same length    
end


times = ones(nData,numel(obj.dataFiles));
eFs = numel(obj.dataFiles)*obj.Fs;   % Equivalent sample rate

% read the instrument frequency into a 2D array
for i = 1: numel(obj.dataFiles)
    %freqs(:,i) = freqs(:,i).*cell2mat(data{1,i}(1:nData)); 
    times(:,i) = times(:,i).*cell2mat(tData{1,i}(1:nData));
end


% We need to align the time since some EPS runs begin slightly ahead or behind the UTC second
tVec = times-floor(times);  % Get a matrix of the fraction of a second, thiese will be tringle waves with 1 minute periods
[r,c] = size(tVec);
lags = zeros(c-1,2*r-1);  % pre-allocate
corr = lags';

% use cross correlation for find the lags
for i = 1:size(tVec,2)-1
    [corr(:,i),lags(i,:)] = xcorr(tVec(:,1),tVec(:,i+1));
end
[~,corrIdx]=max(abs(corr));
maxLags = lags(1,corrIdx);
maxLags = [0,maxLags];
shift = max(maxLags)-maxLags;  % these should show how far each column needs to shift up to time align

% now we are going to go through each of the data file, retreiving and shifting the Fe and RFe
nData = nData-max(shift);
fData = zeros(nData,numel(obj.dataFiles(1)));
rfData = fData;
feCol = find(hdr=='FE');
rfeCol = find(hdr=='RFE');
for i = 1:numel(obj.dataFiles)
    C = readcell(cell2mat(obj.dataFiles(i)));    
    fData(:,i) = cell2mat(C(2+shift(i):nData+shift(i)+1,feCol));
    rfData(:,i) = cell2mat(C(2+shift(i):nData+shift(i)+1,rfeCol));    
end

% Plot the frequency
[t,Fe] = obj.interleaveData(fData,1/obj.Fs,'pos');   % interleave the ETS data
figure(obj.fig); obj.fig=obj.fig+1;
tSettle=settlingTimePlot(t,Fe,1);
ylabel('Frequency error (Hz)')
title(sprintf('Step Frequency Response, Step Response Times = %0.4f s, %0.4f s,',tSettle(1), tSettle(2)))
set(gca,'FontSize',12)

% Plot the ROCOF
[t,RFe] = obj.interleaveData(rfData,1/obj.Fs,'pos');   % interleave the ETS data
figure(obj.fig); obj.fig=obj.fig+1;
tSettle=settlingTimePlot(t,RFe,4);
ylabel('ROCOF error (Hz/s)')
title(sprintf(' Step ROCOF Response, Step Response Times = %0.4f s, %0.4f s,',tSettle(1),tSettle(2)))
set(gca,'FontSize',12)

end


function tSettle = settlingTimePlot(T,Y,effRes)

% the plots step, then step back.  Divide each into two plots and determine the two step respose times
y = [Y(1:floor(length(Y)/2)),Y(floor(length(Y)/2)+1:end)];
t = [T(1:floor(length(T)/2)),T(floor(length(T)/2)+1:end)];
t0 = [0,0];tSettle = [0,0];
for i = 1:2
    yVec = y(:,i);
    tVec = t(:,i);
    [N,edges,~] = histcounts(yVec,100);  % histogrqm
    [~,idx] = max(N);    % index of the highest values in the histogram
    LogicalIndex = yVec>edges(idx) & yVec<=edges(idx+2);
    yMean = mean(yVec(LogicalIndex));
    absYoffset = abs(yMean-yVec);
    %idxT0 = find(absYoffset>effRes,1,'first')-1; % use the effective resolution to determine when the step response begins
    
    % determine a threshold for the step response beginning and end
    thresh = (edges(idx+1)-edges(idx))*effRes;
    idxT0 = find(absYoffset>thresh,1,'first')-1;
    idxSettle = find(absYoffset>thresh,1,'last');
    
    %     % sometimes effRes is not enough so we will multiply it until it is
%     idxSettle = [];
%     thresh = effRes;
%     while (isempty(idxSettle))
%         idxSettle = find(absYoffset>thresh,1,'last');
%         thresh = thresh+effRes;
%     end
    
    
    t0(i) = tVec(idxT0);
    tSettle(i) = tVec(idxSettle);
end

plot(T,Y,'k')
xline(t0(1),'--r');xline(t0(2),'--r');
xline(tSettle(1),'--r');xline(tSettle(2),'--r');
xlabel('Time(s)')

%plot limits
lowLim = t0(1)-0.2; if lowLim<0,lowLim=0;end
hiLim = tSettle(2)+0.2; if hiLim>T(end),hiLim=T(end);end    
xlim([lowLim,hiLim])

tSettle = tSettle-t0;

end

        
