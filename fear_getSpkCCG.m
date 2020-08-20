function fear_getSpkCCG(basename,varargin)

load([basename '.BasicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.Basename '.okUnit.spikes.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
if exist([basicMetaData.Basename '.ripples.events.mat'],'file')
    load([basicMetaData.Basename '.ripples.events.mat'])
else
    ripples.timestamps=[];;
end
load([basicMetaData.Basename '.SleepState.states.mat'])

beh=relabel_ma2sleep(SleepState.MECE.timestamps);

%%

param.binSize=[1,5,25,125]/1000; %sec
param.halfWindow=[0.2,1,5,25]; %sec
param.smoothSD=[1,5,25,125]/1000; %sec

param=parseParameters(param,varargin);

%%
cEdge=unique(okUnit.cluster);
cEdge=[cEdge;max(cEdge)+1];

for tIdx=1:length(param.binSize)
    tBinSize=param.binSize(tIdx);
    fprintf('%s start %s binSize=%d ms\n',datestr(now),basicMetaData.SessionName,tBinSize*1000)
    
    nHalfWin=ceil(param.halfWindow(tIdx)/tBinSize);
    nSmoothWin=ceil(param.smoothSD(tIdx)*4/tBinSize);
    smT=0:nSmoothWin;
    smT=[-fliplr(smT),smT(2:end)]*tBinSize;
    smCore=normpdf(smT,0,param.smoothSD(tIdx));
    smCore=smCore/sum(smCore);
    
    for targetType=1:2
        if targetType==1
            behName='chamber';
            targetName='timestamps';
        else
            behName='homecage';
            targetName='homecage';
        end
        
        for sIdx=1:size(sessions.(targetName),1);
            fprintf('    %s %s-%d\n',datestr(now),behName,sIdx)
            tRange=sessions.(targetName)(sIdx,:);
            
            subClu=okUnit.cluster(okUnit.spikeTime>tRange(1) & okUnit.spikeTime<tRange(2));
            subRes=okUnit.spikeTime(okUnit.spikeTime>tRange(1) & okUnit.spikeTime<tRange(2));
            subBeh=beh(beh(:,2)>tRange(1)&beh(:,2)<tRange(2),:);
            
            clear subState
            if targetType==2
                subState.nrem = subBeh(subBeh(:,3)==3,1:2);
                subState.rem = subBeh(subBeh(:,3)==5,1:2);
                subState.wake = subBeh(subBeh(:,3)==1,1:2);
                
                if ~isempty(ripples.timestamps)
                    subState.swrNrem=ripples.timestamps(...
                        any(ripples.peaks.negative>subState.nrem(:,1)' &...
                        ripples.peaks.negative<subState.nrem(:,2)',2),:);

                    subState.swrWake=ripples.timestamps(...
                        any(ripples.peaks.negative>subState.wake(:,1)' &...
                        ripples.peaks.negative<subState.wake(:,2)',2),:);
                end
            else
                subState=struct();
            end
            subStateList=fieldnames(subState);
            
            for ssIdx=0:length(subStateList)
                if ssIdx==0
                    sName='entire';
                else
                    sName=subStateList{ssIdx};
                end
              
                fprintf('        %s %s\n',datestr(now),sName)

                if ssIdx==0
                    targetRes=subRes;
                    targetClu=subClu;
                else
                    targetRes=subRes(any(subRes>subState.(sName)(:,1)' & subRes<subState.(sName)(:,2)',2));
                    targetClu=subClu(any(subRes>subState.(sName)(:,1)' & subRes<subState.(sName)(:,2)',2));
                end
                
                nSpk=histcounts(targetClu,cEdge);
                
                cnt=CCG([targetRes;-10000],[targetClu;cEdge(end-1)],tBinSize,nHalfWin+nSmoothWin,1);                
                
                smoothed=zeros(nHalfWin*2+1,size(cnt,2),size(cnt,3));
                for cIdx=1:size(cnt,3)
                    temp=Filter0(smCore,cnt(:,:,cIdx));
                    smoothed(:,:,cIdx)=temp(nSmoothWin+1:end-nSmoothWin,:);
                end
                cnt([1:nSmoothWin,end-nSmoothWin+1:end],:,:)=[];
                
                spkCCG.(behName)(sIdx).(sName)(tIdx).raw=cnt;
                spkCCG.(behName)(sIdx).(sName)(tIdx).smoothed=smoothed;
                spkCCG.(behName)(sIdx).(sName)(tIdx).t=(-nHalfWin:nHalfWin)*tBinSize;
                spkCCG.(behName)(sIdx).(sName)(tIdx).nSpk=nSpk;
            end
        end
    end
end
%%
spkCCG.param=param;
spkCCG.generator=mfilename;
spkCCG.generatedate=datestr(now,'yyyy-mm-dd');

%%
save([basicMetaData.AnalysesName '-spkCCG.mat'],'spkCCG','-v7.3')





