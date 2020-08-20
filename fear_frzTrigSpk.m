function  fear_frzTrigSpk(basename,varargin)
% frzTrigSpkHist = fear_frzTrigSpkHist(basename,...)
%  make freeze onset/offset triggered histograms of spikes
%
% options and defaults
%  binSize=0.1; % in sec
%  halfWindow=10; %in sec
%  smoothingSigma=0.5; %in sec
%  minNcell=3;
%  minIsoFrz=10;
%  varName='frzTrigSpkHist';
%  figFileName=[basicMetaData.AnalysesPdf '-frzTrigSpkHist.pdf'];
%  saveFileName=fullfile(basicMetaData.BaseDir,'analyses',[basicMetaData.SessionName '-frzTrigSpkHist.mat']);

%% load basic meta data
load([basename '.basicMetaData.mat']);
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
%% set default
param.tBinSize=0.1; %in sec
param.tBinRange=[-10,10]; %in sec
param.smSD=0.5; %in sec

param.minIsoFrz=5;
%% set options
param=parseParameters(param,varargin);


%% load files

load([basicMetaData.Basename '.okUnit.spikes.mat'])
load([basicMetaData.Basename '.freezeHMM.events.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])
%%

[sesTime,sesNameList]=getChamberTime(basename);
frzDur=diff(freezeHMM.timestamps,1,2);
allFrz.onset=freezeHMM.timestamps(:,1);
allFrz.offset=freezeHMM.timestamps(:,2);
allFrz.isoOnset=allFrz.onset(~any(allFrz.onset>allFrz.offset' &  allFrz.onset<allFrz.offset'+param.minIsoFrz,2));
allFrz.isoOffset=allFrz.offset(~any(allFrz.offset>allFrz.onset'-param.minIsoFrz & allFrz.offset<allFrz.onset',2));
trigList=fieldnames(allFrz);

trigName.onset='all freeze onsets';
trigName.offset='all freeze offsets';
trigName.isoOnset='isolated freeze onsets';
trigName.isoOffset='isolated freeze offsets';

for sesIdx=1:size(sesTime,1)
    for trigIdx=1:length(trigList)
        trigType=trigList{trigIdx};
        frz(sesIdx).(trigType)=allFrz.(trigType)(...
            allFrz.(trigType)>sesTime(sesIdx,1)&...
            allFrz.(trigType)<sesTime(sesIdx,2));
    end
end

%%

smSD=param.smSD;

tBinSize=param.tBinSize;
tBinRange=param.tBinRange;
tBin=tBinRange(1)-4*smSD-tBinSize/2:tBinSize:tBinRange(2)+4*smSD+tBinSize/2;

tBinCenter=(tBin(1:end-1)+tBin(2:end))/2;

cellBin=sort(unique(okUnit.cluster));
cellBin=[cellBin'-0.5,max(cellBin)+0.5];

smBin=0:tBinSize:smSD*4+tBinSize/2;
smBin=[-fliplr(smBin),smBin(2:end)];
smCore=normpdf(smBin,0,smSD);
smCore=smCore/sum(smCore);

for sIdx=1:size(sesTime,1)
    for trigIdx=1:length(trigList)
         trigType=trigList{trigIdx};

         target=frz(sIdx).(trigType);
    
        cnt=zeros([length(tBin),length(cellBin)]-1);
        for tIdx=1:length(target)        
            cnt=cnt+histcounts2(okUnit.spikeTime,okUnit.cluster,tBin+target(tIdx),cellBin);
        end
        fr=cnt/tBinSize/length(target);
        smFR=Filter0(smCore,fr);
    

    
        frzTrigSpk.(trigType).FR{sIdx}=fr(tBinCenter>=tBinRange(1) & tBinCenter<=tBinRange(2),:);
        frzTrigSpk.(trigType).smFR{sIdx}=smFR(tBinCenter>=tBinRange(1) & tBinCenter<=tBinRange(2),:);
        frzTrigSpk.(trigType).trigger{sIdx}=target;
    end
end
frzTrigSpk.chamberName=sesNameList;
frzTrigSpk.generator=mfilename;
frzTrigSpk.generatedate=datestr(now,'yyyy-mm-dd');
frzTrigSpk.param=param;
frzTrigSpk.tBin=tBinCenter(tBinCenter>=tBinRange(1) & tBinCenter<=tBinRange(2)) ;

%% save results

save([basicMetaData.AnalysesName,'-frzTrigSpk.mat'],'frzTrigSpk','-v7.3')



