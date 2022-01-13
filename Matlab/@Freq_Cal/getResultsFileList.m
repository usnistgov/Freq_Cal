function obj = getResultsFileList(obj,varargin)
    %getResultsFileList Prompts the user for the results file path
    % Prompts the user for the path to the Frequency Instrument's results files then gets  lists of results and parameter files

    if isempty(varargin)
        prompt = sprintf('Choose the folder of all results from F0 = %d', obj.F0);
    else
        prompt = sprintf(varargin{1});
    end
    rawDataPath = uigetdir(obj.resultPath,prompt);
    dataNames = obj.getfn(rawDataPath,'.csv');
    paramNames = {};

    % loop through all the data and separate the data from the parameter files
    i = 1;
    while i <= numel(dataNames)
        if contains(dataNames(i),'Parameters')
            paramNames{end+1}=dataNames(i);
            dataNames(i) = [];
        else
            i = i+1;
        end
    end
    obj.paramFiles = paramNames;
    obj.dataFiles = dataNames;
end


% %%-------------------------------------------------------------------------
% % Local functions
% function filenames = getfn(folder,pattern)
%     getfnrec(folder,pattern)
% 
%     idx = ~cellfun(@isempty, regexp(filenames,pattern));
%     filenames =filenames(idx);
% 
%     % This nested function recursively goes through all subfolders
%     % and collects all filenames within them
%         function getfnrec(path,pattern)
%             d = dir(path);
%             filenames = {d(~[d.isdir]).name};
%             filenames = strcat(path,filesep,filenames);
% 
%             dirnames = {d([d.isdir]).name};
%             dirnames = setdiff(dirnames,{'.','..'});
%             for i = 1:numel(dirnames)
%                 fulldirname = [path filesep dirnames{i}];
%                 filenames = [filenames, self.getfn(fulldirname,pattern)];
%             end
%         end
% end