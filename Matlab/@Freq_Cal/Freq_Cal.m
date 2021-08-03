classdef Freq_Cal < handle
    %Top level class file for frequency instrument calibration
    %   
    
    properties
        resultPath
        dataFiles
        paramFiles
        F0 = 50     % Nominal Frequency
        Fs = 50     % reporting rate
        
        fig = 0;    % keep track of open figured
        
        
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

            
            % constructor input arguments
            changed = false;
            for i = 1:2:nargin
                switch varargin{i}
                    case 'F0'
                        obj.F0 = varargin{i+1};
                        changed = true;
                    case 'Fs'
                        obj.Fs = varargin{i+1};
                        changed = true;
%                     case 'Reset'    % delete the .ini file
%                         b = varargin{i+1};
%                         if b=="t" || b=="T"|| b=="true" || b=="True"
%                             if exist(fullfile(appDataPath,name),'file')
%                                 delete(fullfile(appDataPath,name))
%                             end
%                         end
                    otherwise
                        warning('Unrecognized parameter %s',varargin{i});
                end
            end

            if changed
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
      obj = getResultsFileList(obj) % Gets the data and the parameter file list
      calcEffResolution(obj,idx)  % calculates effective resolution from the indexed data and parameter file
      calcFreqRng(obj)  % calcuates the accuracy of frequency and ROCOF for a frequency range test
  end
%%-------------------------------------------------------------------------
  % Private Methods called from external method files
  methods (Access = private)
      
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
        
        
    end
    
end

