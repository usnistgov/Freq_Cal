function calcUnbalancedMagnitude(obj)
% calculated the absolute error from an unbalanced input magnitude test

% Get the measuring range data files
prompt = sprintf('Choose the folder containing the unbalanced magnitude test data');
obj.getResultsFileList(prompt)

% pre-allocate the array using the data from the first file (assume all files are the same number of rows)
C = readcell(cell2mat(obj.dataFiles(1)));
hdr = string(C(1,:));
fCol = find(hdr=='FE');
rfCol = find(hdr=='RFE');
fData = cell2mat(C(2:end,fCol));
rfData = cell2mat(C(2:end,rfCol));
t = linspace(0,length(fData)-1,length(fData))'*1/obj.Fs;

% Plot the frequency error
figure(obj.fig); obj.fig=obj.fig+1;
makePlot(t,fData,obj.MaxAbsFreqError);
ylabel('Frequency absolute error (Hz)')
title('Unbalanced Signal Magnitude Freq Absolute Error')
set(gca,'FontSize',12)

% Plot the ROCOF error
figure(obj.fig); obj.fig=obj.fig+1;
makePlot(t,rfData,obj.MaxAbsRocofError);
ylabel('ROCOF absolute error (Hz)')
title('Unbalanced Signal Magnitude ROCOF Absolute Error')
set(gca,'FontSize',12)


end

function makePlot(t,data,maxErr)
plot(t,abs(data),'k')
yline(maxErr,'--r')
xlabel('Time (s)')
ylim([-0.1*maxErr,(1.1)*maxErr]);
end
