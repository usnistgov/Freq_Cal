classdef Freq_Cal < handle
    %Top level class file for frequency instrument calibration
    %   
    
    properties
        resultPath
        dataFiles
        paramFiles
        F0 = 50     % Nominal Frequency
        Fs = 50     % reporting rate
        
        
        %Delay Times will be calculated
        FreqDelay
        ROCOFDelay
        
        % OpRng and MeasRng are 2 value arrays of doubles:  [lower, upper]
        OpRng
        MeasRng
        
        MaxAbsFreqError
        MaxAbsRocofError
        
        MaxAbsFreqErrorDyn
        MaxAbsRocofErrorDyn
       
        fig = 1;    % keep track of open figures
        
        
    end
    
    %% --------------------------------------------------------------------
    % Constructor
    methods
        function obj = Freq_Cal(varargin)
            %Freq_Cal Construct an instance of this class
            %   obj = Freq_Cal(<option list>)
            %       options: Name,Value pair arguments in any order
            %           "F0" - nominal frequency (default = 50)
            %           "Fs" - instrument reporting rate (measurements per second)
            
            % For user convenience, this will create an ini folder in the user's
            % applications data folder to remember the paths to data and other
            % user parameters
            appDataPath = fullfile(getenv('APPDATA'),'Freq_Cal');
            if ~exist(appDataPath, 'dir' )
                mkdir(appDataPath)
            end
            name = 'Freq_Cal.ini';

            iniExists = exist(fullfile(appDataPath,name),'file');
            if ~iniExists
                obj.resultPath = uigetdir(fullfile(getenv('USERPROFILE'),'Documents'),'Path to instrument data');
            else
                structure = obj.ini2struct(fullfile(appDataPath,name));
                obj.resultPath = structure.ResultsPath.ResultsPath;
                obj.F0 = structure.F0.F0;
                obj.Fs = structure.Fs.Fs;
            end
            
            
            % Optional input parameters
            defaultF0 = 50;
            defaultFs = 50;
            defaultOpRng = [25,75];
            defaultMeasRng = [45,55];
            defaultMaxAbsFreqError = 0.005;
            defaultMaxAbsRocofErr = 0.04;
            defaultMaxAbsFreqErrorDyn = 0.35;
            defaultMaxAbsRocofErrDyn = 14;
            
            p = inputParser;
            
             validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);  
             validRng = @(x) isvector(x) && length(x) == 2 && isa(x,'double');
             
             addParameter(p,'F0',defaultF0,validScalarPosNum)
             addParameter(p,'Fs',defaultF0,validScalarPosNum)
             addParameter(p,'OpRng',defaultOpRng,validRng)
             addParameter(p,'MeasRng',defaultMeasRng,validRng)
             addParameter(p,'MaxAbsFreqError',defaultMaxAbsFreqError,validRng)
             addParameter(p,'MaxAbsRocofError',defaultMaxAbsRocofErr,validRng)
             addParameter(p,'MaxAbsFreqErrorDyn',defaultMaxAbsFreqErrorDyn,validRng)
             addParameter(p,'MaxAbsRocofErrorDyn',defaultMaxAbsRocofErrDyn,validRng)
             
             
             parse(p,varargin{:})
                          
             obj.F0 = p.Results.F0;
             obj.Fs = p.Results.Fs;
             obj.OpRng = p.Results.OpRng;
             obj.MeasRng = p.Results.MeasRng;
             obj.MaxAbsFreqError = p.Results.MaxAbsFreqError;
             obj.MaxAbsRocofError = p.Results.MaxAbsRocofError;           
             obj.MaxAbsFreqErrorDyn = p.Results.MaxAbsFreqErrorDyn;
             obj.MaxAbsRocofErrorDyn = p.Results.MaxAbsRocofErrorDyn;           

            if (obj.F0 ~= defaultF0 && obj.Fs ~= defaultFs)
                structure = struct('ResultsPath',struct('ResultsPath',obj.resultPath),....
                                    'F0',struct('F0',num2str(obj.F0)),...
                                    'Fs',struct('Fs',num2str(obj.Fs)));
                obj.struct2Ini(fullfile(appDataPath,name),structure); % write the .ini file
            end                         
        end
    end

%%-------------------------------------------------------------------------
  % Public Methods called from external method files
  methods (Access = public)
      %obj = getResultsFileList(obj) % Gets the data and the parameter file list
      obj = calcDelayTime(obj);  % calculates the delay time by cross correlation between timestamped reference and instrument readings  
      calcEffResolution(obj,idx)  % calculates effective resolution from the indexed data and parameter file
      calcFreqRng(obj)  % calcuates the accuracy of frequency and ROCOF for a frequency range test
      calcInterHarm(obj) % Calculates accuracy of Frequency and ROCOF under Interharmic tests
      calcDynMeasRange(obj) % Calculates the accuracy during frequency modulation over the measuring range
      calcFreqSettlingTime(obj) % Calculate the settling time of the frequency step using equivalent time sampling
  end
%%-------------------------------------------------------------------------
  % Private Methods called from external method files
  methods (Access = private)
      
  end

  %%  -----------------------------------------------------------------------
    % Local public methods
    methods (Access = public)
        
        function [Freqs, ROCOFs] = shiftByDelayTime(obj,Freqs,ROCOFs)  
            % shifts series of frequency and ROCOF measurements by the estimated delay times
            
        end
        
        
    end
  
%%  -----------------------------------------------------------------------
    % Local Static Methods
    methods(Static)
        
        % write a structure to a .ini file
        % Dirk Lohse (2021). struct2ini (
        % https://www.mathworks.com/matlabcentral/fileexchange/22079-struct2ini), 
        % MATLAB Central File Exchange. Retrieved January 15, 2021.
        function struct2Ini(file,structure)
            fid = fopen(file,'w');
            
            sects = fieldnames(structure);
            
            for i = 1:numel(sects)
                sect = char(sects(i));
                fprintf(fid,'\n[%s]\n',sect);
                mem = structure.(sect);
                if ~isempty(mem)
                    memNames = fieldnames(mem);
                    for j = 1:numel(memNames)
                        memName = char(memNames(j));
                        memVal = structure.(sect).(memName);
                        fprintf(fid,'%s=%s\n',memName,memVal);
                    end                    
                end
            end
            fclose(fid);
        end
        
        % read a structure from a .ini file
        %freeb (2021). ini2struct
        %(https://www.mathworks.com/matlabcentral/fileexchange/45725-ini2struct),
        %MATLAB Central File Exchange. Retrieved January 15, 2021.
        function Struct = ini2struct(FileName)
            % Parses .ini file
            % Returns a structure with section names and keys as fields.
            %
            % Based on init2struct.m by Andriy Nych
            % 2014/02/01
            f = fopen(FileName,'r');                    % open file
            while ~feof(f)                              % and read until it ends
                s = strtrim(fgetl(f));                  % remove leading/trailing spaces
                if isempty(s) || s(1)==';' || s(1)=='#' % skip empty & comments lines
                    continue
                end
                if s(1)=='['                            % section header
                    Section = genvarname(strtok(s(2:end), ']'));
                    Struct.(Section) = [];              % create field
                    continue
                end
                
                [Key,Val] = strtok(s, '=');             % Key = Value ; comment
                Val = strtrim(Val(2:end));              % remove spaces after =
                
                if isempty(Val) || Val(1)==';' || Val(1)=='#' % empty entry
                    Val = [];
                elseif Val(1)=='"'                      % double-quoted string
                    Val = strtok(Val, '"');
                elseif Val(1)==''''                     % single-quoted string
                    Val = strtok(Val, '''');
                else
                    Val = strtok(Val, ';');             % remove inline comment
                    Val = strtok(Val, '#');             % remove inline comment
                    Val = strtrim(Val);                 % remove spaces before comment
                    
                    [val, status] = str2num(Val);       
                    if status, Val = val; end           % convert string to number(s)
                end
                
                if ~exist('Section', 'var')             % No section found before
                    Struct.(genvarname(Key)) = Val;
                else                                    % Section found before, fill it
                    Struct.(Section).(genvarname(Key)) = Val;
                end
            end
            fclose(f);                        
        end
        
        
        function filenames = getfn(folder,pattern)
            getfnrec(folder,pattern)
            
            idx = ~cellfun(@isempty, regexp(filenames,pattern));
            filenames =filenames(idx);
            
            % This nested function recursively goes through all subfolders
            % and collects all filenames within them
            function getfnrec(path,pattern)
                d = dir(path);
                filenames = {d(~[d.isdir]).name};
                filenames = strcat(path,filesep,filenames);
                
                dirnames = {d([d.isdir]).name};
                dirnames = setdiff(dirnames,{'.','..'});
                for i = 1:numel(dirnames)
                    fulldirname = [path filesep dirnames{i}];
                    filenames = [filenames, self.getfn(fulldirname,pattern)];
                end
            end
        end
        
        function [t, Y] = interleaveData(Data,dT,dir)
            % Interleave 2D ETS data in a positive (first column first) or
            % negative (last column first) direction
            
            % time vector
            [nSamples,nColumns] = size(Data);
            t_axis = (0:1:nSamples-1)*dT;
            t_tot = zeros(nSamples,nColumns);
            for i = 1:nColumns
               t_tot(:,i) = t_axis' + ((i-1)*dT/nColumns); 
            end
            t_final = t_tot;
            if dir == 'pos', t_final = fliplr(t_final); end
            t_final = t_final(:);
            [t,idx] = sort(t_final);
            
            % Data vector
            Y_tot = Data;
            %if dir == 'pos',Y_tot = fliplr(Y_tot); end
            Y_final = Y_tot(:);
            Y = Y_final(idx);                                    
        end
        
    end
    
end

