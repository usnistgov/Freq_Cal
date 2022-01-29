function calcRocofSettlingTime(obj)
% calculates the settling time using instrument data taken during a frequency step
% The data is from a series of tests with the step time decremented by 1/10 reporting cycle
% This function reads the instrumnt data, interleaves the data and calculated the 10%
% and 90% boundaries and the settling time.

% Get the measuring range data files
prompt = sprintf('Choose the folder containing the ROCOF Settling Time data');
obj.getResultsFileList(prompt)

% pre-allocate the array using the data from the first file (assume all files are the same number of rows)
C = readcell(cell2mat(obj.dataFiles(1)));
hdr = string(C(1,:));
col = find(hdr=='Inst ROCOF');

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


rocofs = ones(nData,numel(obj.dataFiles));
eFs = numel(obj.dataFiles)*obj.Fs;   % Equivalent sample rate

% read the instrument frequency into a 2D array
for i = 1: numel(obj.dataFiles)
    rocofs(:,i) = rocofs(:,i).*cell2mat(data{1,i}(1:nData));        
end
[t,Y] = obj.interleaveData(rocofs,1/obj.Fs,'pos');   % interleave the ETS data


% Use a Histogram to find the levels
[N,edges,bin] = histcounts(Y,100);  % histogrqm
[~,iFirst] = max(N);    % index of the highest values in the histogram
N(iFirst) = 0;
[~,iSecond] = max(N);   % index of the second highest values in the histogram

% now get the mean values of the beginning and ending states
firstLogicalIndex = Y>edges(iFirst) & Y<=edges(iFirst+1);
meanFirst = mean(Y(firstLogicalIndex));

secondLogicalIndex = Y>edges(iSecond) & Y<=edges(iSecond+1);
meanSecond = mean(Y(secondLogicalIndex));

%next we need to determine if the step is positive or negative
stepSize = abs(meanFirst-meanSecond);
halfStep = min([meanFirst,meanSecond])+(stepSize/2);
idx = find(Y>halfStep,1,'first'); % find the index of the first high sample
posStep = false;
if idx > length(Y)/4,posStep = true;end % postive or negative Step

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
secondHighLimit = meanSecond + .05*stepSize;
secondLowLimit = meanSecond - .05*stepSize;

midLevel = (meanFirst+meanSecond)/2; % The value at the middle of the step


%find the step
idxMid = find(abs(Y)>=abs(midLevel),1,'first');
t = t-t(idxMid);  % This is the delay corrected time

idxStart = idxMid - (25/obj.F0 * eFs);
idxEnd = idxMid + (25/obj.F0 * eFs);


% idxFirst = find(firstLogicalIndex,1,'last');
% idxSecond = find(secondLogicalIndex,1,'first');
% idxStart = idxFirst - (14/obj.F0 * eFs);
% idxEnd = idxSecond + (14/obj.F0 * eFs);
% 
% if posStep
%     idxMid = find(Y>=midLevel,1,'first');
% else
%     idxMid = find(Y<=midLevel,1,'first');
% end
% 
% t = t-t(idxMid);  % This is the delay corrected time

% trim the data and the time vector
Y = Y(idxStart:idxEnd);
t = t(idxStart:idxEnd);

figure(obj.fig);obj.fig = obj.fig+1;
plot(t,Y,'*')




 if posStep
     tSecond = t(Y > secondLowLimit);
 else    
     tSecond = t(Y < secondHighLimit); 
 end

 % % ARG:  At the time of this writing, I have not created a fitting routine
% to determine the refence values of a frequency step, so I will use the
% effective resolution to determine when the step response begins.  this
% may actually produce a slightly shorter settling time since some of the
% delay time (such as ethernet polling intervals) may not be included in the reference.
% This may be better becaues polling intervals can be non deterministic

% find the first Y value greater than the meanFirst +/- Effective Resolution
idxPlus = find(Y>=meanFirst+obj.RocofEffRes,1,'first');
idxMinus = find(Y<=meanFirst-obj.RocofEffRes,1,'first');
idx = min([idxPlus,idxMinus])/eFs+t(1);
refFreq = (t < idx)*meanFirst + (t>=idx)*meanSecond;
firstST1 = idx;

% create limit lines
% line = ones(length(tFirst),1);
% firstHighLine = line*firstHighLimit;
% firstLowLine = line*firstLowLimit;

line = ones(length(tSecond),1);
secondHighLine = line*secondHighLimit;
secondLowLine = line*secondLowLimit;


% find the last place thet Y crosses the second high or low limit lines
tLastLow = t(find(Y<=secondLowLimit,1,'last'));
tLastHigh = t(find(Y>=secondHighLimit,1,'last'));
lastST1 = max([tLastHigh,tLastLow]);
ST1 = lastST1 - firstST1;


resid = abs(Y - refFreq);
idx = find(resid>obj.MaxAbsFreqError,1,'first');
firstST2 = t(idx);
idx = find(resid>obj.MaxAbsFreqError,1,'last');
lastST2 = t(idx);

ST2 = lastST2 - firstST2;


% plot
figure(obj.fig); obj.fig = obj.fig+1;
%subplot(2,1,1)
% p =plot(...
%      t,refFreq,'g',...
%      t,Y,'k',...
%      tFirst,firstHighLine,'r--',...
%      tFirst,firstLowLine,'r--',...
%      tSecond,secondHighLine,'r--',...
%      tSecond,secondLowLine,'r--');

% we need to "uncorrect" the time step (remove the delay time)
if posStep
    iStep = find(refFreq>=midLevel,1,'first');
else
    iStep = find(refFreq<=midLevel,1,'first');
end
tStep = t(iStep);

    


p =plot(...
     t-tStep,refFreq,'g',...
     t-tStep,Y,'k',...
     tSecond-tStep,secondHighLine,'r--',...
     tSecond-tStep,secondLowLine,'r--');


 set(p(1),'color',[0,0.5,0]);
 
 %xline(firstST1,'b')
 xline(lastST1-tStep,'b')

 % legend position
 if posStep
     legPos = 'southeast';
 else
     legPos = 'northeast';
 end
 legend('Reference','Instrument','location',legPos)
 
 xlabel('Time from frequency step of input (s)');
 ylabel('Frequency (Hz)');
 title(sprintf('Settling time using 95%% of the step = %0.4f s',ST1))
 axis padded
 xlim([t(1),t(end)])
 set(gca,'FontSize',12)
 
%  subplot(2,1,2)
%  plot(t,abs(resid),'k')
%  set(gca,'YScale', 'log')
%  
%  yline(obj.MaxAbsFreqError,'r--')
%  xline(firstST2,'b')
%  xline(lastST2,'b')
% 
%  xlabel('Time from step (s)');
%  ylabel('Log absolute value of error (Hz)')
%  title(sprintf('Settling time using absolute error = %0.4f s',ST2))
% 
%  xlim([t(1),t(end)])
%  ylim([-1,6])
end

