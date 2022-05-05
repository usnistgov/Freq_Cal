function calcDropAndRecover(obj)
% calculates the settling time using instrument data taken during a frequency step
% The data is from a series of tests with the step time decremented by 1/10 reporting cycle
% This function reads the instrumnt data, interleaves the data and calculated the 10%
% and 90% boundaries and the settling time.

% Get the measuring range data files
prompt = sprintf('Choose the folder containing the Noise Test data');
obj.getResultsFileList(prompt)

% pre-allocate the array using the data from the first file (assume all files are the same number of rows)
C = readcell(cell2mat(obj.dataFiles(1)));
hdr = string(C(1,:));

fCol = find(hdr=='FE');
rfCol = find(hdr=='RFE');
fData = cell2mat(C(2:end,fCol));
rfData = cell2mat(C(2:end,rfCol));
t = linspace(0,length(fData)-1,length(fData))'*1/obj.Fs;

% Plot the frequency
figure(obj.fig); obj.fig=obj.fig+1;
settlingTimePlot(t,fData,obj.MaxAbsFreqError);
ylabel('Absolute Frequency Error (Hz)')
title('Noise Test Abs Frequency Error')
set(gca,'FontSize',12)

% Plot the ROCOF
figure(obj.fig); obj.fig=obj.fig+1;
settlingTimePlot(t,rfData,obj.MaxAbsRocofError);
ylabel('Absolute ROCOF Error (Hz/s)')
title('Noise Test ABS ROCOF Error')
set(gca,'FontSize',12)

end


function settlingTimePlot(t,Y,errMax)
plot(t,abs(Y),'k')
yline(errMax,'--r')
xlabel('Time(s)')
xlim([0,5])
ylim([0,errMax+.001])













end

        
% N(iFirst) = 0;
% [~,iSecond] = max(N);   % index of the second highest values in the histogram
% 
% % now get the mean values of the beginning and ending states
% firstLogicalIndex = Y>edges(iFirst) & Y<=edges(iFirst+1);
% meanFirst = mean(Y(firstLogicalIndex));
% 
% secondLogicalIndex = Y>edges(iSecond) & Y<=edges(iSecond+1);
% meanSecond = mean(Y(secondLogicalIndex));
% 
% %next we need to determine if the step is positive or negative
% stepSize = abs(meanFirst-meanSecond);
% halfStep = min([meanFirst,meanSecond])+(stepSize/2);
% 
% 
% idx = find(Y>halfStep,1,'first'); % find the index of the first high sample
% posStep = false;
% if idx > length(Y)/4,posStep = true;end % postive or negative Step
% 
% % next find out if the most common value matches the initial state or not
% if posStep && (meanFirst < meanSecond)
%     ;
% elseif ~posStep && (meanSecond < meanFirst)
%     ;
% else  % swap the values if it does not match
%     temp = firstLogicalIndex;
%     firstLogicalIndex = secondLogicalIndex;
%     secondLogicalIndex = temp;
%     temp = meanFirst;
%     meanFirst = meanSecond;
%     meanSecond=temp;
% end
% 
% stepSize = abs(meanFirst - meanSecond);
% firstHighLimit = meanFirst + .05*stepSize;
% firstLowLimit = meanFirst - .05*stepSize;
% secondHighLimit = meanSecond + .05*stepSize;
% secondLowLimit = meanSecond - .05*stepSize;
% 
% midLevel = (meanFirst+meanSecond)/2; % The value at the middle of the step
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
% 
% t = t-t(idxMid);  % This is the delay corrected time
% 
% % trim the data and the time vector
% Y = Y(idxStart:idxEnd);
% t = t(idxStart:idxEnd);
% 
%  if posStep
% %     tFirst = t(Y < firstHighLimit);
%      tSecond = t(Y > secondLowLimit);
%  else    
% %     tFirst = t(Y > firstLowLimit);
%      tSecond = t(Y < secondHighLimit); 
%  end
% 
% % % Settling time using 95% of the step
% % firstST1 = tFirst(end);
% % lastST1 = tSecond(1);
% % ST1 = lastST1 - firstST1;
% 
% % % ARG:  At the time of this writing, I have not created a fitting routine
% % to determine the refence values of a frequency step, so I will use the
% % effective resolution to determine when the step response begins.  this
% % may actually produce a slightly shorter settling time sinse some of the
% % delay time may not be included in the reference
% 
% % find the first Y value greater than the meanFirst +/- Effective Resolution
% idxPlus = find(Y>meanFirst+obj.FreqEffRes,1,'first');
% idxMinus = find(Y<meanFirst-obj.FreqEffRes,1,'first');
% idx = min([idxPlus,idxMinus])/eFs+t(1);
% refFreq = (t < idx)*meanFirst + (t>=idx)*meanSecond;
% firstST1 = idx;
% 
% % create limit lines
% % line = ones(length(tFirst),1);
% % firstHighLine = line*firstHighLimit;
% % firstLowLine = line*firstLowLimit;
% 
% line = ones(length(tSecond),1);
% secondHighLine = line*secondHighLimit;
% secondLowLine = line*secondLowLimit;
% 
% 
% % find the last place thet Y crosses the second high or low limit lines
% tLastLow = t(find(Y<=secondLowLimit,1,'last'));
% tLastHigh = t(find(Y>=secondHighLimit,1,'last'));
% lastST1 = max([tLastHigh,tLastLow]);
% ST1 = lastST1 - firstST1;
% 
% 
% resid = abs(Y - refFreq);
% idx = find(resid>obj.MaxAbsFreqError,1,'first');
% firstST2 = t(idx);
% idx = find(resid>obj.MaxAbsFreqError,1,'last');
% lastST2 = t(idx);
% 
% ST2 = lastST2 - firstST2;
% 
% 
% % plot
% figure(obj.fig); obj.fig = obj.fig+1;
% %subplot(2,1,1)
% % p =plot(...
% %      t,refFreq,'g',...
% %      t,Y,'k',...
% %      tFirst,firstHighLine,'r--',...
% %      tFirst,firstLowLine,'r--',...
% %      tSecond,secondHighLine,'r--',...
% %      tSecond,secondLowLine,'r--');
% 
% % we need to "uncorrect" the time step (remove the delay time)
% if posStep
%     iStep = find(refFreq>=midLevel,1,'first');
% else
%     iStep = find(refFreq<=midLevel,1,'first');
% end
% tStep = t(iStep);
% 
%     
% 
% 
% p =plot(...
%      t-tStep,refFreq,'g',...
%      t-tStep,Y,'k',...
%      tSecond-tStep,secondHighLine,'r--',...
%      tSecond-tStep,secondLowLine,'r--');
% 
% 
%  set(p(1),'color',[0,0.5,0]);
%  
%  %xline(firstST1,'b')
%  xline(lastST1-tStep,'b')
% 
%  % legend position
%  if posStep
%      legPos = 'southeast';
%  else
%      legPos = 'northeast';
%  end
%  legend('Reference','Instrument','location',legPos)
%  
%  xlabel('Time from frequency step of input (s)');
%  ylabel('Frequency (Hz)');
%  title(sprintf('Settling time using 95%% of the step = %0.4f s',ST1))
%  axis padded
%  xlim([t(1),t(end)])
%  set(gca,'FontSize',12)
%  
% %  subplot(2,1,2)
% %  plot(t,abs(resid),'k')
% %  set(gca,'YScale', 'log')
% %  
% %  yline(obj.MaxAbsFreqError,'r--')
% %  xline(firstST2,'b')
% %  xline(lastST2,'b')
% % 
% %  xlabel('Time from step (s)');
% %  ylabel('Log absolute value of error (Hz)')
% %  title(sprintf('Settling time using absolute error = %0.4f s',ST2))
% % 
% %  xlim([t(1),t(end)])
% %  ylim([-1,6])
% end
