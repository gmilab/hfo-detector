classdef hfodat
    properties
        nsamp    % num samples
        nchan    % num channels
        fs       % sampling rate
        data     % data as [samples] x [channels]
        time     % vector of time points (length of nsamp)
        label    % cell array of channel labels
        
        trigs    % cell array [nchan] x 1 each containing vector of gold standard / known HFO locations in samples
        
        markers  % optional: more detailed info
        
        dataset_name
    end
    
    methods
        function pstr = output_ripplelab(obj, idx)
            pstr.data = obj.data(:,idx);
            pstr.fs = obj.fs;
        end
        
        function hfodata = run_ripplelab(obj, ch_idx, method)
            % run the specified ripplelab algorithm on the specified channel(s)
            
            addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'RIPPLELAB', 'Functions')));
            addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'RIPPLELAB', 'External', 'wafo')));
            addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'RIPPLELAB', 'External', 'freetb4matlab')));
            addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'RIPPLELAB', 'External', 'Misc')));
            %       addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'RIPPLELAB', 'Functions', 'HFO'));
            addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'external', 'ripplelab'));
            
            params = struct();
            params.s_FreqIni = 80;
            params.s_FreqEnd = 512;
            
            switch method
                case 'hil'
                    params.s_EpochTime = 60;
                    params.s_SDThres = 5;
                    params.s_MinWind = 10;
                    
                    fn = @rl_findHFOxHIL;
                    
                case 'mni'
                    params.s_EpochTime = 60;
                    params.s_EpoCHF = 60;
                    params.s_PerCHF = 95;
                    params.s_MinWin = 10;
                    params.s_MinGap = 10;
                    params.s_ThresPerc = 99.9999;
                    params.s_BaseSeg = 125;
                    params.s_BaseShift = 0.5000;
                    params.s_BaseThr = 0.6700;
                    params.s_BaseMin = 5;
                    
                    fn = @rl_findHFOxMNI;
                    
                case 'sll'
                    params.s_EpochTime = 60;
                    params.s_FiltWind = 5;
                    params.s_Percentil = 97.5000;
                    params.s_MinWind = 12;
                    params.s_FiltEq = 0;
                    
                    fn = @rl_findHFOxSLL;
                    
                case 'ste'
                    params.s_EpochTime = 60;
                    params.s_RMSWindow = 3;
                    params.s_RMSThres = 5;
                    params.s_MinWind = 6;
                    params.s_MinTime = 10;
                    params.s_NumOscMin = 6;
                    params.s_BPThresh = 3;
                    
                    fn = @rl_findHFOxSTE;
                otherwise
                    error('Unknown method');
            end
            
            hfodata = struct();
            hfodata.trigs = cell(length(ch_idx), 1);
            for kk = 1:length(ch_idx)
                hfodata.trigs{kk} = fn(obj.data(:,ch_idx(kk)), obj.fs, params);
                hfodata.trigs{kk} = mean(hfodata.trigs{kk}, 2);
            end
            hfodata.channels = ch_idx;
            hfodata.detector = method;
            hfodata.label = obj.label(ch_idx);
            hfodata.dataset_name = obj.dataset_name;
        end
        
        function hfodata = run_zurich(obj, varargin)
            % RUN_ZURICH - Run the Zurich detector on
            % hfodata = dat.run_zurich()
            % hfodata = dat.run_zurich(channels)
            
            ip = inputParser;
            addOptional(ip, 'ch_idx', 1:obj.nchan);
            addParameter(ip, 'parameter_set', 'FRSpecPara', @(x) any(strcmp(x, {'FRSpecPara', 'RSpecPara'})));
            parse(ip, varargin{:});
            
            ZurichDetectorPath = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'ZurichDetector');
            addpath(ZurichDetectorPath)
            
            hfo = Core.HFO;
            paras = load(fullfile(ZurichDetectorPath, '+Demo', 'Spec', 'ECoG', 'Parameters', ...
                sprintf('%s.mat', ip.Results.parameter_set)));
            hfo.Para = paras.DetPara;
            hfo.Data = [];
            hfo.Data.signal = obj.data(:,ip.Results.ch_idx);
            hfo.Data.channelNames = obj.label(ip.Results.ch_idx);
            hfo.Data.sampFreq = obj.fs;
            hfo.Data.maxIntervalToJoin = hfo.Para.maxIntervalToJoinPARA * obj.fs;
            hfo.Data.MinHighEntrIntvLen = hfo.Para.MinHighEntrIntvLenPARA * obj.fs;
            hfo.Data.minEventTime = hfo.Para.minEventTimePARA * obj.fs;
            hfo.Data.sigDurTime = diff(obj.time([1 end]));
            hfo.Data.timeInterval = [1 120];  % baseline interval
            hfo.Data.nbChannels = obj.nchan;
            hfo.Data.nbSamples = obj.nsamp;
            hfo.Data.dataSetup = [];
            
            smoothBool = false;
            hfo = getFilteredSignal(hfo, smoothBool);
            
            hfo = getBaseline(hfo);
            
            RefType   = 'spec';
            hfo = getEvents(hfo, RefType);
            disp(hfo.Events)
            
            hfo = getEventProperties(hfo);
            disp(hfo.Events)
            
            hfo = getRefinements(hfo, RefType);
            disp(hfo.Refinement)
            
            maskPassMultiChan  = hfo.Refinement.maskMultChanRefine;
            CondMask = hfo.Refinement.maskEventCondSelect;
            Events = hfo.Events;
            hfo = hfo.RefineEvents(Events,1,CondMask, maskPassMultiChan);
            
            hfodata = struct();
            
            hfodata.trigs = cell(length(ip.Results.ch_idx), 1);
            for kk = 1:length(ip.Results.ch_idx)
                hfodata.trigs{kk} = hfo.RefinedEvents{1}.Markings.start{kk} + (hfo.RefinedEvents{1}.Markings.len{kk} ./ 2);
            end
            hfodata.channels = ip.Results.ch_idx;
            hfodata.detector = 'ZurichNCH';
            hfodata.label = obj.label(ip.Results.ch_idx);
            hfodata.dataset_name = obj.dataset_name;
        end
        
        function hfodata = run_hscdetector(obj, varargin)
            addpath /d/gmi/1/simeon/hfo_detector/repo/MATLAB
            addpath /d/gmi/1/simeon/hfo_detector/repo/MATLAB/fieldtrip
            
            ip = inputParser;
            addParameter(ip, 'channels', 'all');
            parse(ip, varargin{:});
            
            ch_idx = ip.Results.channels;
            if strcmp(ch_idx, 'all')
                ch_idx = 1:obj.nchan;
            end
            
            % convert to fieldtrip format
            data_cnt = [];
            data_cnt.trial = {obj.data'};
            data_cnt.fsample = obj.fs;
            data_cnt.time = {0:(1/obj.fs):((obj.nsamp-1)/obj.fs)};
            data_cnt.label = obj.label;
            
            hfodata = struct();
            hfodata.trigs = cell(length(ch_idx), 1);
            for kk = 1:length(ch_idx)
                hfodata.trigs{kk} = ieeg_hfodetector(data_cnt, 'channel', ch_idx(kk));
            end
            hfodata.channels = ch_idx;
            hfodata.detector = 'hsc';
            hfodata.label = obj.label(ch_idx);
            hfodata.dataset_name = obj.dataset_name;
        end
        
        function confusion_data = compute_confusion(obj, hfodata, varargin)
            % take output from run_ripplelab and compare to the known gold standard in the object
            ip = inputParser;
            addParameter(ip, 't_threshold', 0.25);
            parse(ip, varargin{:});
            
            % initialize table
            confusion_data = table();
            confusion_data.label = obj.label(hfodata.channels);
            confusion_data.detector = repmat({hfodata.detector}, length(hfodata.channels), 1);
            
            % initialize descriptives
            confusion_data.TP = zeros(height(confusion_data),1);
            confusion_data.FP = zeros(height(confusion_data),1);
            confusion_data.FN = zeros(height(confusion_data),1);
            confusion_data.PPV = zeros(height(confusion_data),1);
            confusion_data.SEN = zeros(height(confusion_data),1);
            
            % loop over channels and compute PPV and SEN
            for kk = 1:length(hfodata.channels)
                ch = hfodata.channels(kk);
                
                gold_standard = obj.trigs{ch};
                thishfotrigs = hfodata.trigs{kk};
                
                for nn = 1:length(thishfotrigs)
                    if any(abs(gold_standard - thishfotrigs(nn)) <= (obj.fs * ip.Results.t_threshold))
                        % confusion_data.TP(kk) = confusion_data.TP(kk) + 1;
                        % >>> edit: count true positives on the gold standard instead
                    else
                        confusion_data.FP(kk) = confusion_data.FP(kk) + 1;
                    end
                end
                
                for nn = 1:length(gold_standard)
                    if any(abs(gold_standard(nn) - thishfotrigs) <= (obj.fs * ip.Results.t_threshold))
                        confusion_data.TP(kk) = confusion_data.TP(kk) + 1;
                    else
                        confusion_data.FN(kk) = confusion_data.FN(kk) + 1;
                    end
                end
                
                sum_TP = confusion_data.TP(kk);
                sum_FP = confusion_data.FP(kk);
                sum_FN = confusion_data.FN(kk);
                
                confusion_data.PPV(kk) = sum_TP/(sum_TP+sum_FP);
                confusion_data.SEN(kk) = sum_TP/(sum_TP+sum_FN);
            end
        end
        
        function confusion_data = run_all_ripplelab(obj, varargin)
            ip = inputParser;
            addParameter(ip, 'detector', {'sll', 'ste', 'hil', 'mni'});
            addParameter(ip, 'channels', 1:obj.nchan);
            parse(ip, varargin{:});
            
            confusion_data = [];
            for kk = 1:length(ip.Results.detector)
                hfodata = obj.run_ripplelab(ip.Results.channels, ip.Results.detector{kk});
                ct = obj.compute_confusion(hfodata);
                
                if isempty(confusion_data)
                    confusion_data = ct;
                else
                    confusion_data = cat(1, confusion_data, ct);
                end
            end
        end
    end
    
    methods(Static)
        function summary = summarize_confusion(confusion_data)
            if any(contains(confusion_data.label, '0dB'))
                labels = {' 0dB'; ' 5dB'; ' 10dB'; ' 15dB'};
            elseif any(contains(confusion_data.label, 'snr'))
                labels = {'snr9'; 'snr6'; 'snr3'; 'snr0'};
            end
            detectors = unique(confusion_data.detector);
            
            % generate mesh
            [labels, detectors] = meshgrid(labels, detectors);
            
            % initialize table
            summary = table();
            summary.label = reshape(labels, [], 1);
            summary.detector = reshape(detectors, [], 1);
            
            summary.TP = zeros(height(summary),1);
            summary.FP = zeros(height(summary),1);
            summary.FN = zeros(height(summary),1);
            summary.PPV = zeros(height(summary),1);
            summary.SEN = zeros(height(summary),1);
            
            % loop through conditions, compute sum
            for kk = 1:height(summary)
                rows_idx = contains(confusion_data.label, summary.label{kk}) & strcmp(confusion_data.detector, summary.detector{kk});
                
                summary.TP(kk) = sum(confusion_data.TP(rows_idx), 1);
                summary.FP(kk) = sum(confusion_data.FP(rows_idx), 1);
                summary.FN(kk) = sum(confusion_data.FN(rows_idx), 1);
                
                summary.PPV(kk) = summary.TP(kk)/(summary.TP(kk)+summary.FP(kk));
                summary.SEN(kk) = summary.TP(kk)/(summary.TP(kk)+summary.FN(kk));
            end
        end
        
        function obj = load_donos()
            obj = hfodat();
            obj.data = [];
            
            for kk = 1:2
                % dsnum is either 1 or 2
                load(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data', 'donos',...
                    sprintf('simulation_data%d.mat', kk)), 'databip');
                
                % label with dsnum
                labels = databip.label(1:4);
                labels = cellfun(@(x) sprintf('%d-%s', kk, x), labels, 'UniformOutput', false);
                
                obj.data = cat(2, obj.data, databip.channel(:,1:4));
                obj.label = cat(1, obj.label, labels);
            end
            
            obj.nsamp = size(databip.channel,1);
            obj.nchan = size(obj.data,2);
            
            obj.fs = databip.fs;
            obj.time = 0:(1/databip.fs):((obj.nsamp-1)/databip.fs);
            
            % HFOs are at the same place in all channels
            obj.trigs = repmat({1000:2000:obj.nsamp}, obj.nchan, 1);
            obj.dataset_name = 'donos';
        end
        
        function obj = load_anywave()
            obj = hfodat();
            
            % there's only one anywave dataset
            filename = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data', 'anywave', 'simulation_rate_3.ades');
            
            [filedir,fileroot,~] = fileparts(filename);
            
            if ~exist(fullfile(filedir, [fileroot, '.ades']), 'file')
                error('%s does not exist!', fullfile(filedir, [fileroot, '.ades']))
            end
            
            if ~exist(fullfile(filedir, [fileroot, '.dat']), 'file')
                error('%s does not exist!', fullfile(filedir, [fileroot, '.dat']))
            end
            
            %% parse header
            fid = fopen(fullfile(filedir, [fileroot, '.ades']), 'r');
            hdr_text = textscan(fid, '%s', 'Delimiter', '\n', 'HeaderLines', 1);
            fclose(fid);
            
            hdr_text = hdr_text{1};
            hdr_text = cellfun(@(x) strtrim(strsplit(x, '=')), hdr_text, 'UniformOutput', false);
            hdr_text = cat(1, hdr_text{:});
            
            % get sampling rate and nsamples
            idx = strcmp(hdr_text(:,1), 'samplingRate');
            obj.fs = str2double(hdr_text{idx, 2});
            
            idx = strcmp(hdr_text(:,1), 'numberOfSamples');
            obj.nsamp = str2double(hdr_text{idx, 2});
            
            hdr_text = hdr_text(3:end,:);
            obj.nchan = size(hdr_text,1);
            
            obj.label = hdr_text(:,1);
            
            %% load data
            fid = fopen(fullfile(filedir, [fileroot, '.dat']), 'rb');
            obj.data = fread(fid, 'float32');
            fclose(fid);
            
            obj.data = reshape(obj.data, obj.nchan, obj.nsamp);
            obj.data = obj.data';
            
            obj.time = 0:(1/obj.fs):((obj.nsamp-1)/obj.fs);
            
            %% parse marker file
            if exist(fullfile(filedir, [fileroot, '.mrk']), 'file')
                obj.markers = readtable(fullfile(filedir, [fileroot, '.mrk']), 'FileType', 'text', 'Delimiter', '\t');
                obj.markers.Properties.VariableNames(1:4) = {'Label', 'Value', 'Time', 'Duration'};
                
                if length(obj.markers.Properties.VariableNames) == 5
                    obj.markers.Properties.VariableNames{5} = 'Misc';
                end
            else
                obj.markers = [];
            end
            obj.markers.HFO = contains(obj.markers.Label, 'Ripple');
            
            %% construct trigs cell
            obj.trigs = cell(obj.nchan, 1);
            for kk = 1:obj.nchan
                obj.trigs{kk} = obj.markers.Time(obj.markers.HFO & contains(obj.markers.Misc, [obj.label{kk},',']));
                obj.trigs{kk} = round(obj.trigs{kk} .* obj.fs);    % convert to samples
            end
            
            obj.dataset_name = sprintf('anywave');
        end
    end
end