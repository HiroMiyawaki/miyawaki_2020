function varargout = fear_frzTrigSpkHist_KS2(basename,varargin)

load([basename '.basicMetaData.mat']);
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
%% set default
param.binSize=[0.01,0.1,1]; %in sec
param.halfWindow=[1,10,100]; %in sec
param.smoothingSigma=[0.05,0.5,5]; %in sec
param.minNcell=3;
param.minIsoFrz=10;
param.varName='frzTrigSpkHist';
param.figFileName=[basicMetaData.AnalysesPdf '-frzTrigSpkHist.pdf'];
param.saveFileName=fullfile(basicMetaData.BaseDir,'analyses',[basicMetaData.SessionName '-frzTrigSpkHist.mat']);
%% set options
param=parseParameters(param,varargin);


%% load files

load([basicMetaData.Basename '.okUnit.spikes.mat'])
load([basicMetaData.Basename '.freezeHMM.events.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])
%%

sesTime=sessions.timestamps;
% split Cue and Extinction in thirds
extIdx=find(strcmpi({basicMetaData.chamber.name},'CueAndExtinction'));
FirstPip=min(cues.timestamps.Pip(cues.timestamps.Pip(:,1)>=sesTime(extIdx,1),1));

temp=cues.timestamps.Tone(cues.timestamps.Tone(:,1)>=sesTime(extIdx,1),1);
extPipStart=min(cues.timestamps.Pip(cues.timestamps.Pip(:,1)>=temp(9)));

borders=[sesTime(extIdx,1),FirstPip,extPipStart,sesTime(extIdx,2)];

sesTime=[sesTime(1:extIdx-1,:)
    borders(1:end-1)',borders(2:end)'
    sesTime(extIdx+1:end,:)];

sesNameList={basicMetaData.chamber(1:extIdx-1).name};

sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Base'];
sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Cue'];
sesNameList{end+1}=[basicMetaData.chamber(extIdx).name '-Ext'];
sesNameList={sesNameList{:},basicMetaData.chamber(extIdx+1:end).name};


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



%% get cell list to be used
okCells=1:length(okUnit.cluInfo.channel);
%% get region list to be plotsn
okUnit.cluInfo.region(contains(okUnit.cluInfo.region,{'PrL L2','PrL L3','PrL L1'}))={'PrL L2/3'};
regList=unique(okUnit.cluInfo.region);
nCell=sum(strcmp(repmat(okUnit.cluInfo.region,length(regList),1),repmat(regList',1,length(okCells))),2)';

regList=regList(nCell>=param.minNcell);

%% get bin numbers and smoothing core
for binIdx=1:length(param.binSize)
    nSmoothWindow(binIdx)=ceil(param.smoothingSigma(binIdx)*3/param.binSize(binIdx));
    
    smoothingCore{binIdx}=normpdf((-nSmoothWindow(binIdx):nSmoothWindow(binIdx))*param.binSize(binIdx),0,param.smoothingSigma(binIdx));
    smoothingCore{binIdx}=smoothingCore{binIdx}/sum(smoothingCore{binIdx});
    
    nHalfBin(binIdx)=ceil(param.halfWindow(binIdx)/param.binSize(binIdx));
end
%% get PETH
clear frzTrigSpkHist
for regIdx=1:length(regList)
    regName=strrep(regList{regIdx},' ','_');
    regName=strrep(regName,'/','_');
    target=find(strcmp(okUnit.cluInfo.region,regList{regIdx}));
    
    fprintf('%s processing %d cells in %s (%d/%d regions)\n',datestr(now),length(target),regList{regIdx},regIdx,length(regList))
    
    for trigIdx=1:length(trigList)
        for binIdx=1:length(param.binSize)
            trigType=trigList{trigIdx};
            frzTrigSpkHist.(regName).(trigType)(binIdx).fr=zeros(length(target),2*nHalfBin(binIdx)+1,length(frz));
            frzTrigSpkHist.(regName).(trigType)(binIdx).smoothFr=zeros(length(target),2*nHalfBin(binIdx)+1,length(frz));
            frzTrigSpkHist.(regName).(trigType)(binIdx).z=zeros(length(target),2*nHalfBin(binIdx)+1,length(frz));
            frzTrigSpkHist.(regName).(trigType)(binIdx).smoothZ=zeros(length(target),2*nHalfBin(binIdx)+1,length(frz));
            frzTrigSpkHist.(regName).(trigType)(binIdx).nTrigger(sesIdx)=length(frz(sesIdx).(trigType));
        end
    end
    
    cellCnt=0;
    for cIdx=target
        cellCnt=cellCnt+1;
        spkTime=okUnit.spikeTime(okUnit.cluster==cIdx);
        
        for binIdx=1:length(param.binSize)
            fr=hist(spkTime,basicMetaData.detectionintervals.lfp(1):param.binSize(binIdx):basicMetaData.detectionintervals.lfp(2))/param.binSize(binIdx);
            fr(1)=[];fr(end)=[];
            meanFR(binIdx)=mean(fr);
            stdFR(binIdx)=std(fr);
        end
        for sesIdx=1:length(frz);
            for binIdx=1:length(param.binSize)
                subset=spkTime(spkTime>sesTime(sesIdx,1)-param.halfWindow(binIdx) & spkTime<sesTime(sesIdx,2)+param.halfWindow(binIdx));
                for trigIdx=1:length(trigList)
                    trigType=trigList{trigIdx};
                    if isempty(frz(sesIdx).(trigType))
                        fr=nan(2*nHalfBin(binIdx)+1,1);
                        smoothed=nan(2*nHalfBin(binIdx)+1,1);
                    else
                        [cnt,t]=CCG([subset;frz(sesIdx).(trigType)],...
                            [ones(size(subset));2*ones(size(frz(sesIdx).(trigType)))],...
                            param.binSize(binIdx),nHalfBin(binIdx)+nSmoothWindow(binIdx),1);
                        
                        fr=squeeze(cnt(nSmoothWindow(binIdx)+1:end-nSmoothWindow(binIdx),1,end))/length(frz(sesIdx).(trigType))/param.binSize(binIdx);
                        smoothed=conv(squeeze(cnt(:,1,end))/length(frz(sesIdx).(trigType))/param.binSize(binIdx),smoothingCore{binIdx},'same');
                        smoothed=smoothed(nSmoothWindow(binIdx)+1:end-nSmoothWindow(binIdx));
                    end
                    
                    
                    frzTrigSpkHist.(regName).(trigType)(binIdx).fr(cellCnt,:,sesIdx)=fr;
                    frzTrigSpkHist.(regName).(trigType)(binIdx).z(cellCnt,:,sesIdx)=(fr-meanFR(binIdx))/stdFR(binIdx);
                    frzTrigSpkHist.(regName).(trigType)(binIdx).smoothFr(cellCnt,:,sesIdx)=smoothed;
                    frzTrigSpkHist.(regName).(trigType)(binIdx).smoothZ(cellCnt,:,sesIdx)=(smoothed-meanFR(binIdx))/stdFR(binIdx);
                    frzTrigSpkHist.(regName).(trigType)(binIdx).t=t(nSmoothWindow+1:end-nSmoothWindow)/1e3;
                end
            end
        end
    end
    
end
%
%% make summary pdf
if ~isempty(param.figFileName)
    close all
    pageCnt=0;
    doAppend='';
    outputList={};
    tempDir='~/Desktop/tempForFrzTrigSpkHist';
    
    if exist(tempDir,'dir')
        cleanTempDir=false;
    else
        mkdir(tempDir);
        cleanTempDir=true;
    end
    
    sortSessionName='CueAndExtinction-Cue';
    for binIdx=1:length(param.binSize)
    for trigIdx=1:length(trigList)
        trigType=trigList{trigIdx};
        
        fh=initFig4A4('landscape',true,'fontsize',7);
        pageCnt=pageCnt+1;
        for regIdx=1:length(regList);
            regName=strrep(regList{regIdx},' ','_');
            regName=strrep(regName,'/','_');
            temp=frzTrigSpkHist.(regName).(trigType)(binIdx);
            temp.smoothFr(isnan(temp.smoothFr))=0;
            
            sortIdx=find(strcmpi(sesNameList,sortSessionName));
            
            [maxRate,peakPos]=max(squeeze(temp.smoothFr(:,:,sortIdx)),[],2)    ;
            peakPos(maxRate==0)=inf;
            [~,order]=sort(peakPos);
            
            for sesIdx=1:size(temp.smoothFr,3)
                subplot2(length(regList)+1,size(temp.smoothFr,3),regIdx,sesIdx)
                scale=max(squeeze(temp.smoothFr(order,:,sesIdx)),[],2);
                scale(scale==0)=1;
                imagesc(temp.t,1:size(temp.smoothFr,1), squeeze(temp.smoothFr(order,:,sesIdx))./scale)
                box off
                if sesIdx==1;ylabel(regList{regIdx});end
                if regIdx==1;title(sesNameList{sesIdx});end
                if regIdx==length(regList);
                    if contains(trigType,'onset','IgnoreCase',true)
                        xlabel(['Time from freeze onset(s)']);
                    else
                        xlabel(['Time from freeze offset(s)']);
                    end
                end
                colormap(gca,hot)
                ax=fixAxis;
            end
        end
        addScriptName(mfilename)
        subplot2(length(regList)+1,size(temp.smoothFr,3),length(regList)+1,1)
        axis off
        ax=fixAxis;
        text2(-0.5,0.75,[basicMetaData.SessionName ' ' trigName.(trigType) ' triggered spike histogram - bin size =' num2str(param.binSize(binIdx)) ' s'],...
            ax,'fontsize',12,'verticalAlign','top')
        
        txt={};
        switch lower(trigType)
            case 'onset'
                txt{1}='All freeze onsets are used as triggers';
            case 'offset'
                txt{1}='All freeze offsets are used as triggers';
            case 'isoonset'
                txt{1}=sprintf('Only freeze onsets isolated > %d sec from previous Freeze offsets are used as triggers',param.minIsoFrz);
            case 'isooffset'
                txt{1}=sprintf('Only freeze offsets isolated > %d sec from next Freeze onsets are used as triggers',param.minIsoFrz);
        end
        
        txt{2}='Only OK cells are used';
        
        txt{3}=sprintf('Counted in %0.1f sec bin and smoothed with Gaussian (\\sigma = %0.1f sec)',param.binSize,param.smoothingSigma);
        txt{4}='Cue & Extinction session are separated into baseline(before first tone),Cue (1st-8th tone),and Ext (9th tone-)';
        txt{5}='Cells are sorted by peak position of Cue&Extinction-Cue';
        
        text2(0,0.5,txt,ax,'verticalAlign','top')
        
        outputList{end+1}=fullfile(tempDir,['page' num2str(pageCnt)]);
                print(fh,outputList{end},'-dpdf','-painters','-r300')
    end
    end
    
        catPDF(param.figFileName,outputList,'remove',true)
    
    if cleanTempDir
        rmdir(tempDir,'s')
    end
end

%% save results
for regIdx=1:length(regList)
    regName=strrep(regList{regIdx},' ','_');
    regName=strrep(regName,'/','_');
    frzTrigSpkHist.(regName).param=param;
    frzTrigSpkHist.(regName).generator=mfilename;
    frzTrigSpkHist.(regName).generatedate=today('datetime');
end
if ~isempty(param.varName) &&  ~isempty(param.saveFileName)
    if strcmp(param.varName,'frzTrigSpkHist')
        eval([param.varName '=frzTrigSpkHist;']);
    end
    save(param.saveFileName,param.varName,'-v7.3')
end

%% set output
if nargout >0
    varargout{1}=frzTrigSpkHist;
end

