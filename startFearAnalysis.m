clear

rootDir='~/data/Fear/triple';

% set first bit '1' to include the session into the analyses
% set second bit '1' to redo analyses
% set third bit '1' to redo sleep scoring

sessionList={
     1,0,0,'achel180320'
     1,0,0,'booyah180430'
     1,0,0,'chimay180612'
     1,0,0,'duvel190505'
     1,0,0,'estrella180808'
     1,0,0,'feuillien180830'
     1,0,0,'guiness181002'
     1,0,0,'hoegaarden181115'
     1,0,0,'innis190601'
     1,0,0,'jever190814'
     1,0,0,'karmeliet190901'
     1,0,0,'leffe200124'
     1,0,0,'maredsous200224'
     1,0,0,'nostrum200304'
     1,0,0,'oberon200325'
             };


         

%%
clc
for sessionListIndex=1:size(sessionList,1)
    if sessionList{sessionListIndex,1}~=1; continue; end
    
    clearvars -except rootDir sessionList sessionListIndex       
    sessionName=sessionList{sessionListIndex,end}; %session name
    reAnalyses=sessionList{sessionListIndex,2}==1;
    reScoring=sessionList{sessionListIndex,3}==1;
    

    %when _AnalysesMetadataText haven't generated, only detect shock
    if ~exist(fullfile(rootDir,sessionName,[sessionName '_AnalysesMetadataText.m']),'file')
        basicMetaData.Basename=fullfile(rootDir,sessionName,sessionName);
        basicMetaData.dat=fullfile('~/data/Fear/triple_raw/dat/',[sessionName '.dat']);
        basicMetaData.SessionName=sessionName;
        basicMetaData.BaseDir=fullfile(rootDir,sessionName);
        basicMetaData.nCh=206;
        basicMetaData.SampleRates.dat=20e3;
        shockTTLch=202:204;
        info=dir(basicMetaData.dat);
        basicMetaData.nSample.dat=info.bytes/basicMetaData.nCh/2;
        detectionintervals=[0,info.bytes/basicMetaData.nCh/2/basicMetaData.SampleRates.dat];
        detectionintervals=detectionintervals+30*[1,-1]; %ignore first and last 30 sec
        if ~exist([basicMetaData.Basename '.shocks.events.mat'],'file')
            if ~exist([basicMetaData.BaseDir '/lfp/'],'dir')
                fprintf('make %s\n', [basicMetaData.BaseDir '/lfp/'])
                mkdir([basicMetaData.BaseDir '/lfp/']);
            end
            detectTTLpulses(basicMetaData,...
                'chList',shockTTLch,...
                'minInterval',[0.1,1e-4,1e-4],...
                'minDuration',[0.1,1e-4,1e-4],...
                'detectionintervals',detectionintervals,...
                'filetype','dat',...
                'evtFileName',[basicMetaData.BaseDir '/lfp/' basicMetaData.SessionName '.shk.evt'],...
                'saveFileName',[basicMetaData.Basename '.shocks.events.mat'],...
                'varName','shocks',...
                'chName',{ 'ShockTrig'    'ShockL'    'ShockR'});    
        end
        
        continue
    end
    
    
    fprintf('\n\n%s start %s \n\n',datestr(now),sessionName)
    
    run(fullfile(rootDir,sessionName,[sessionName '_AnalysesMetadataText.m']));    
    
    fearAnalysesPipeline(rootDir, sessionName,chName,detectionintervals,...
                    'ttlCh',ttlCh, 'accelerometerCh',accelerometerCh,...
                    'videofile',videofile, 'chamber',chamber,...
                    'reAnalyses',reAnalyses, 'detectionCh',detectionCh)
    
    %sleep scoring
    if reScoring || ~exist(fullfile(rootDir,sessionName,[sessionName '.SleepState.states.mat']),'file')
        input=prepareForTheStateEditor(rootDir,sessionName); 
        TheStateEditor(fullfile(rootDir,sessionName,sessionName),input)
        uiwait;
        addMECEtoSleepState(basicMetaData.Basename)
    end
    
    %do routine for spikes
    if  (exist(fullfile(rootDir, sessionName,'ks2','probe1'),'dir') && ...
        exist(fullfile(rootDir, sessionName,'ks2','probe2'),'dir') && ...
        exist(fullfile(rootDir, sessionName,'ks2','probe3'),'dir'))|| ...
        exist(fullfile(rootDir, sessionName,[sessionName '.okUnit.spikes.mat']),'file')
        
       
        % analyses for spikes
        fearAnalysesSpikePipeline(rootDir, sessionName,'reAnalyses',reAnalyses)
        
        % analyses for reactivated ensembles
        fearAnalysesReactPipeline(rootDir, sessionName,'reAnalyses',reAnalyses)
        
        fprintf('\n\n%s %s is done\n\n\n',datestr(now),sessionName)
    end

end


%% making figures
coactPaper_fig1
coactPaper_fig2
coactPaper_fig3
coactPaper_fig4
coactPaper_fig5

coactPaper_figS01
coactPaper_figS02
coactPaper_figS03
coactPaper_figS04
coactPaper_figS05
coactPaper_figS06
coactPaper_figS07
coactPaper_figS08
coactPaper_figS09
coactPaper_figS10
coactPaper_figS11
coactPaper_figS12
coactPaper_figS13
coactPaper_figS14
coactPaper_figS15
coactPaper_figS16
coactPaper_figS17
coactPaper_figS18
coactPaper_figS19

coactPaper_tableS01
coactPaper_tableS02
coactPaper_tableS03
coactPaper_tableS04



