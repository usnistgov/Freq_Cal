function calcModulation(obj)
% This is some experimentation with the modulation raw data

% Get the data file
prompt = sprintf('Choose the folder containing the modulation data');
rawDataPath = uigetdir(obj.resultPath,prompt);
dataNames = obj.getfn(rawDataPath,'.mat');
obj.dataFiles = dataNames;

j = 1;
load(obj.dataFiles{j})
[~,Name,~]=fileparts(obj.dataFiles{j});
obj.Fs = P(1).SampleRate;

i = 1;
k = 1;

Y = P(i).Samples(k,:).';
N = length(Y);
tn = (-(N/2):(N/2)-1).'/obj.Fs;

Data = timeseries(Y,tn,'Name',Name);
Data = setuniformtime(Data,'Interval',1/obj.Fs);

HHT = HilbertHuang_class('TimeSeries',Data,'EMD',true, 'Hilbert',true);

HHT.plot('IMFs');
HHT.plot('Hilbert');

end