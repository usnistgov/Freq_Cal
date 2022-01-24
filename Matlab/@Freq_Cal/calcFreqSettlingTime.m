function calcFreqSettlingTime(obj)
% calculates the settling time using instrument data taken during a frequency step
% The data is from a series of tests with the step time decremented by 1/10 reporting cycle
% This function reads the instrumnt data, interleaves the data and calculated the 10%
% and 90% boundaries and the settling time.

% Get the measuring range data files
prompt = sprintf('Choose the folder containing the Freqency Settling Time data');
obj.getResultsFileList(prompt)

% pre-allocate the array using the data from the first file (assume all files are the same number of rows)
C = readcell(cell2mat(obj.dataFiles(1)));
hdr = string(C(1,:));
col = find(hdr=='Inst. Freq');

% read all the data into a cell array.  some data is longer than others
data = cell(1,numel(obj.dataFiles));
C = readcell(cell2mat(obj.dataFiles(1)));
data{1,1} = C(2:end,col);
nData = length(data{1,1});

for i = 2:numel(obj.dataFiles)
    C = readcell(cell2mat(obj.dataFiles(i)));
    data{1,i} = C(2:end,col);
    if length(data{1,i})<nData,nData = length(data{1,i});end
end


freqs = ones(nData,numel(obj.dataFiles));

% read the instrument frequency into a 2D array
for i = 1: numel(obj.dataFiles)
    freqs(:,i) = freqs(:,i).*cell2mat(data{1,i}(1:nData));        
    %freqs(:,i) = freqs(:,i).*cell2mat(data{1,i}((length(data{1,i})-nData+1):end)); 
end
[t,Y] = obj.interleaveData(freqs,1/obj.Fs,'pos');   % interleave the ETS data

idxMid = find(Y>obj.F0,1,'first'); % find the index of the first high sample
posStep = false;
if idxMid > length(Y)/4,posStep = true;end % postive or negative Step

% use a Histogram to find the 2 most common data values
[N,edges,bin] = histcounts(Y,100);  % histogrqm
[~,iFirst] = max(N);    % index of the highest values in the histogram
N(iFirst) = 0;
[~,iSecond] = max(N);   % index of the second highest values in the histogram

% now get the mean values of the beginning and ending states
firstLogicalIndex = Y>edges(iFirst) & Y<=edges(iFirst+1);
meanFirst = mean(Y(firstLogicalIndex));

secondLogicalIndex = Y>edges(iSecond) & Y<=edges(iSecond+1);
meanSecond = mean(Y(secondLogicalIndex));

% next find out if the most common value matches the initial state or not
if posStep && (meanFirst < meanSecond)
    ;
elseif ~posStep && (meanSecond < meanFirst)
    ;
else  % swap the values if it does not match
    temp = firstLogicalIndex;
    firstLogicalIndex = secondLogicalIndex;
    secondLogicalIndex = temp;
    temp = meanFirst;
    meanFirst = meanSecond;
    meanSecond=temp;
end

stepSize = abs(meanFirst - meanSecond);
firstHighLimit = meanFirst + .1*stepSize;
firstLowLimit = meanFirst - .1*stepSize;
secondHighLimit = meanSecond + .1*stepSize;
secondLowLimit = meanSecond - .1*stepSize;

midLevel = (meanFirst+meanSecond)/2; % The value at the middle of the step

posStep = meanFirst < meanSecond;
if posStep
    idxFirst = find(firstLogicalIndex,1,'last');
    idxSecond = find(secondLogicalIndex,1,'first');
    idxStart = idxFirst - (14/obj.F0 * 10 * obj.Fs);
    idxEnd = idxSecond + (14/obj.F0 * 10 * obj.Fs);
    idxMid = find(Y>=midLevel,1,'first');
else
    idxFirst = find(firstLogicalIndex,1,'first');
    idxSecond = find(secondLogicalIndex,1,'last');
    idxStart = idxSecond - (14/obj.F0 * 10 * obj.Fs);
    idxEnd = idxFirst + (14/obj.F0 * 10 * obj.Fs);
    idxMid = find(Y<=midLevel,1,'last');     
end

t = t-t(idxMid);

% trim the data and the time vector
Y = Y(idxStart:idxEnd);
t = t(idxStart:idxEnd);

if posStep
    tFirst = t(Y < firstHighLimit);
    tSecond = t(Y > secondLowLimit);
else    
    tFirst = t(Y < firstHighLimit);
    tSecond = t(Y > secondLowLimit); 
end

% Settling time using 90% of the step
firstST1 = tFirst(end);
lastST1 = tSecond(1);
ST1 = lastST1 - firstST1;

% create limit lines
line = ones(length(tFirst),1);
firstHighLine = line*firstHighLimit;
firstLowLine = line*firstLowLimit;

line = ones(length(tSecond),1);
secondHighLine = line*secondHighLimit;
secondLowLine = line*secondLowLimit;

% Create a reference data set
refFreq = (t < 0).*meanFirst + (t>=0).*meanSecond;

resid = Y - refFreq;
idx = find(resid>obj.MaxAbsFreqError,1,'first');
firstST2 = t(idx);
idx = find(resid>obj.MaxAbsFreqError,1,'last');
lastST2 = t(idx);

ST2 = lastST2 - firstST2;


% plot
figure(obj.fig); obj.fig = obj.fig+1;
subplot(2,1,1)
plot(...
     t,refFreq,'g',...
     t,Y,'k',...
     tFirst,firstHighLine,'r--',...
     tFirst,firstLowLine,'r--',...
     tSecond,secondHighLine,'r--',...
     tSecond,secondLowLine,'r--')
 
 xline(firstST1,'b')
 xline(lastST1,'b')

 legend('Reference','Instrument','location','southeast')
 xlabel('Time from step (s)');
 ylabel('Frequency (Hz)');
 title(sprintf('Settling time using 90%% of the step = %0.4f s',ST1))
 xlim([t(1),t(end)])
 
 subplot(2,1,2)
 plot(t,abs(resid),'k')
 set(gca,'YScale', 'log')
 
 yline(obj.MaxAbsFreqError,'r--')
 xline(firstST2,'b')
 xline(lastST2,'b')

 xlabel('Time from step (s)');
 ylabel('Log absolute value of error (Hz)')
 title(sprintf('Settling time using absolute error = %0.4f s',ST2))

 xlim([t(1),t(end)])
 ylim([-1,6])
end
