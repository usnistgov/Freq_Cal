function calcRamp(obj)
% calculated the absolute error from an unbalanced input magnitude test

% Get the measuring range data files
prompt = sprintf('Choose the folder containing the frequency ramp test data');
obj.getResultsFileList(prompt)

% Each pair of result files constitutes one test scenario
i = 1;
while i < numel(obj.dataFiles)
    fData = [];fRefData = [];rfData = [];rfRefData = [];
    settleSamples = ceil(obj.FreqSettle*obj.Fs);  % the number of samples to ignore ant the beginning and end
    first = ceil(2+3*settleSamples);
    
    % frequency data    
    C = readcell(cell2mat(obj.dataFiles(i)));     
    last = length(C)- 3*settleSamples;
    hdr = string(C(1,:));
    fData = cell2mat(C(first:last,find(hdr=='FE')));
    fRefData = cell2mat(C(first:last,find(hdr=='REF_FREQ')));
    
    % ROCOF data
    %settleSamples = ceil(obj.RocofSettle*obj.Fs);  % the number of samples to ignore ant the beginning and end    
    rfData = cell2mat(C(first:last,find(hdr=='RFE')));
    rfRefData = cell2mat(C(first:last,find(hdr=='REF_ROCOF')));
        
    i = i+1;  
    
    % frequency data    
    C = readcell(cell2mat(obj.dataFiles(i)));     
    hdr = string(C(1,:));
    fData = vertcat(fData,cell2mat(C(first:last,find(hdr=='FE'))));
    fRefData = vertcat(fRefData,cell2mat(C(first:last,find(hdr=='REF_FREQ'))));
    
    % ROCOF data
    %settleSamples = ceil(obj.RocofSettle*obj.Fs);  % the number of samples to ignore ant the beginning and end    
    rfData = vertcat(rfData,cell2mat(C(first:last,find(hdr=='RFE'))));
    
    i = i+1;
    
    
    figure(obj.fig); obj.fig=obj.fig+1;
    plot(abs(fData),'k');
    x1 = 0:100:numel(fData)-1;
    set(gca,'XTick',x1,'XTickLabel',round(fRefData(1:100:end)*10)/10);
    ylabel('Absolute frequency error (Hz)')
    xlabel('Frequency (Hz)')
    %yline(obj.MaxAbsFreqError,'--r')
    
    figure(obj.fig); obj.fig=obj.fig+1;
    plot(abs(rfData),'k');
    x1 = 0:100:numel(fData)-1;
    set(gca,'XTick',x1,'XTickLabel',round(fRefData(1:100:end)*10)/10);
    ylabel('Absolute ROCOF error (Hz)')
    xlabel('Frequency (Hz)')
    %yline(obj.MaxAbsRocofError,'--r')
    

end

end
% % pre-allocate the array using the data from the first file (assume all files are the same number of rows)
% C = readcell(cell2mat(obj.dataFiles(1)));
% fCol = find(hdr=='FE');
% rfCol = find(hdr=='RFE');
% fData = cell2mat(C(2:end,fCol));
% rfData = cell2mat(C(2:end,rfCol));
% t = linspace(0,length(fData)-1,length(fData))'*1/obj.Fs;
% 
% % Plot the frequency error
% figure(obj.fig); obj.fig=obj.fig+1;
% makePlot(t,fData,obj.MaxAbsFreqError);
% ylabel('Frequency absolute error (Hz)')
% title('Unbalanced Signal Magnitude Freq Absolute Error')
% set(gca,'FontSize',12)
% 
% % Plot the ROCOF error
% figure(obj.fig); obj.fig=obj.fig+1;
% makePlot(t,rfData,obj.MaxAbsRocofError);
% ylabel('ROCOF absolute error (Hz)')
% title('Unbalanced Signal Magnitude ROCOF Absolute Error')
% set(gca,'FontSize',12)
% 
% 
% end
% 
% function makePlot(t,data,maxErr)
% plot(t,abs(data),'k')
% yline(maxErr,'--r')
% xlabel('Time (s)')
% ylim([-0.1*maxErr,(1.1)*maxErr]);
% end
