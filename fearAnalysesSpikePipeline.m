function fearAnalysesSpikePipeline(rootDir,sessionName,varargin)
load(fullfile(rootDir,sessionName,[sessionName '.basicMetaData.mat']))    

reAnalyses=false;
for idx=1:length(varargin)/2
    varName=varargin{2*idx-1};
    varVal=varargin{2*idx};

    switch lower(varName)
        case lower('reAnalyses')
            reAnalyses=varVal;
        otherwise
            warning(['wrong option name : ' varName])
            return                
    end
end

%% make spike file based on KS2
if ~exist([basicMetaData.AnalysesName '-ks2.spikes.mat'],'file')
    datFile=fullfile('~/Desktop/',basicMetaData.SessionName,[basicMetaData.SessionName '.dat']);
    if ~exist(datFile,'file')
        datFile=basicMetaData.dat;
    end
    makeSpikeFiles_KS2(basicMetaData.Basename,datFile)
end

% get fet and waveform based on KS2
if ~exist([basicMetaData.AnalysesName,'-waveform.mat'],'file')
%     datFile=fullfile('/Volumes/temporal/forKK2/',basicMetaData.SessionName,[basicMetaData.SessionName '.dat']);
    datFile=fullfile('~/Desktop/',basicMetaData.SessionName,[basicMetaData.SessionName '.dat']);
    getFet_KS2(basicMetaData.Basename,datFile)    
end

% get spike stats based on KS2
if ~exist([basicMetaData.AnalysesName '-clusterStats.mat'],'file')
    getClusterQuality_KS2(basicMetaData.Basename)
end

% get spike waveform stats based on KS2
if ~exist([basicMetaData.AnalysesName '-waveformStats.mat'],'file')
    getWaveformStats_KS2(basicMetaData.Basename)
end

%% select ok units based on cluster stats
if  ~exist([basicMetaData.Basename '.okUnit.spikes.mat'],'file')
    getOKunit_KS2(basicMetaData.Basename);    
end

%% separate cluinfo from **unit.spikes
if reAnalyses || ~exist([basicMetaData.AnalysesName '-okUnit.cellinfo.mat'],'file')
    
    fear_makeUnitInfo(basicMetaData.Basename)
    
end

%% get reactivations first
if reAnalyses || ~exist([basicMetaData.AnalysesName '-icaReac.mat'],'file')
    fear_icaReac(basicMetaData.Basename)
end


%% get mean FR of each behavioral state
if reAnalyses || ~exist([basicMetaData.AnalysesName '-meanFR.mat'],'file')
    fear_getMeanFR(basicMetaData.Basename);
end

%% freeze modulation
if reAnalyses || ~exist(fullfile(basicMetaData.BaseDir,'analyses',[basicMetaData.SessionName '-frzFR.mat']),'file')
    fear_FR_freezeModulation_KS2(basicMetaData.Basename,...
        'minWakeDur',40,'minIsoDist',15,'maxIsiIdx',0.1,'minNcell',3,'minFR',1e-4,...
        'showPlot',true,'varName','frzFR',...
        'figFileName',[basicMetaData.AnalysesPdf '-frzFR.pdf'],...
        'saveFileName',[basicMetaData.AnalysesName '-frzFR.mat']);
end

if reAnalyses || ~exist([basicMetaData.AnalysesName '-frzFRmod.mat'],'file')
    fear_getfrzFRmod(basicMetaData.Basename)
end



%% freeze triggered spike histograms
if reAnalyses || ~exist(fullfile(basicMetaData.BaseDir,'analyses',[basicMetaData.SessionName '-frzTrigSpkHist.mat']),'file')
    fear_frzTrigSpkHist_KS2(basicMetaData.Basename,...
        'figFileName',[basicMetaData.AnalysesPdf '-frzTrigSpkHist.pdf'],...
        'saveFileName',[basicMetaData.AnalysesName '-frzTrigSpkHist.mat'],...
        'varName','frzTrigSpkHist',...
        'minNcell',3,'minIsoFrz',10);
end

if reAnalyses || ~exist([basicMetaData.AnalysesName,'-frzTrigSpk.mat'],'file')
    fear_frzTrigSpk(basicMetaData.Basename)
end

%% detect synaptic connections
if reAnalyses || ~exist([basicMetaData.AnalysesName '-synConn.mat'],'file')
    fprintf('%s start %s with data of %s\n',datestr(now),'detectSynapticConnections',basicMetaData.SessionName)
    detectSynapticConnections(basicMetaData.Basename, '.okUnit.spikes.mat','nSurrogate',1000,'spkVar','okUnit')
end

%% SWR triggered histogram
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-swrTrigHist.mat'])) && ...
        exist([basicMetaData.Basename '.ripples.events.mat'],'file')
    fear_getSWR_spike_PETH(basicMetaData.Basename);
end

%subdevide each homecage sessions
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-swrTrigHistSub.mat'])) && ...
        exist([basicMetaData.Basename '.ripples.events.mat'],'file')
        fear_getSwrTrigHist_SubdivideHomecage(basicMetaData.Basename);
end

%subdivide each homecage sessions with longer window
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-swrTrigHistSubLong.mat'])) && ...
        exist([basicMetaData.Basename '.ripples.events.mat'],'file')
        fear_getSwrTrigHist_SubdivideHomecage(basicMetaData.Basename,...
        'tBinSize',10e-3,'nHalfBin',500,'tSmoothSigma',100e-3,...
        'varName','swrTrigHistSubLong','matFileName',[basicMetaData.AnalysesName '-swrTrigHistSubLong.mat']);
end


%% spindle triggered spike histogram
if reAnalyses || ~exist([basicMetaData.AnalysesName '-spindleTrigHist.mat'],'file')
    fear_getSpindle_spike_PETH(basicMetaData.Basename)
end

%% amy HFO triggered histogram
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-hfoTrigHist.mat'])) && ...
        exist([basicMetaData.Basename '.amyHFO.events.mat'],'file')
     fear_getHFO_spike_PETH(basicMetaData.Basename);
end

%% spike CCG
if reAnalyses || ~exist([basicMetaData.AnalysesName '-spkCCG.mat'],'file')
    fear_getSpkCCG(basicMetaData.Basename)
end

if ~exist([basicMetaData.AnalysesPdf ,'-fear_homecage_spkCCG_1ms.pdf'],'file')
    fear_homecage_spkCCG_1ms(basicMetaData.Basename)
end

%% classify excitatory/inihbitory cells
if exist([basicMetaData.AnalysesName '-okUnit.cellinfo.mat'],'file')
    load([basicMetaData.AnalysesName '-okUnit.cellinfo.mat'])
    if ~isfield(okUnitInfo,'cellType')
        if exist([basicMetaData.AnalysesName '-synConn.mat'],'file')
            load([basicMetaData.AnalysesName '-synConn.mat'])
            if ~isfield(synConn,'inspected')
                synnConnEditor(basicMetaData.Basename)
                uiwait;
            end
            fear_EIseparation(basicMetaData.Basename)
        end
    end
end
%% cue triggered spike hist
if reAnalyses || ~exist([basicMetaData.AnalysesName '-cueTrigSpk.mat'],'file')
    fprintf('%s getting cue triggered spike histogram\n',datestr(now))
    fear_cueTrigSpk(basicMetaData.Basename)
end

if reAnalyses || ~exist([basicMetaData.AnalysesName '-tonePETHeach.mat'],'file')
    fprintf('%s getting cue triggered firing rates for each tone/cell \n',datestr(now))
    fear_tonePETHeach(basicMetaData.Basename)
end
%%
if (reAnalyses || ~exist([basicMetaData.AnalysesName '-ripFrMod.mat'],'file'))&& ...
         exist([basicMetaData.Basename '.ripples.events.mat'],'file')
    fear_getRipFrMod(basicMetaData.Basename)
end
if reAnalyses || ~exist([basicMetaData.AnalysesName '-spdlFrMod.mat'],'file')
    fear_getSpdlFrMod(basicMetaData.Basename)
end
if reAnalyses || ~exist([basicMetaData.AnalysesName '-deltaFrMod.mat'],'file')
    fear_getDeltaFrMod(basicMetaData.Basename)
end
if reAnalyses || ~exist([basicMetaData.AnalysesName '-hfoFrMod.mat'],'file')
    fear_getHfoFrMod(basicMetaData.Basename)
end

