function fear_shockTrigCoact(basename,varargin)
%%
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

%%
param.binSize=0.02; % in sec
param.halfWindow=5; % in sec
param.targetHC=3;

param.varName='evtTrigCoactStr';
param.saveFile=[basicMetaData.AnalysesName '-shockTrigCoact.mat'];

param.reacFile=[basicMetaData.AnalysesName '-icaCoactTimeCondHT.mat'];
%%
param=parseParameters(param,varargin);

%%
load([basicMetaData.Basename '.SleepState.states.mat'])
load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.shocks.events.mat'])

temp=load(param.reacFile);
vName=fieldnames(temp);
coactTime=temp.(vName{1});
%%
temp=[shocks.timestamps.ShockL;shocks.timestamps.ShockR];
for n=1:size(shocks.timestamps.ShockTrig,1)
    idx=find(temp(:,1)>shocks.timestamps.ShockTrig(n,1),1,'first');
    trig(n)=temp(idx);
end

    tRange=sessions.timestamps(2,:);
        
    target=1:length(coactTime.sigLevel);
    
    time=[];
    val=[];
    
    nWin=ceil(param.halfWindow/param.binSize);
    tBin=(-nWin:nWin)*param.binSize;
    tBorder=[-inf,((-nWin:nWin+1)-0.5)*param.binSize,inf];

    each.rate=nan(length(target),length(trig),nWin*2+1);
    each.strength=nan(length(target),length(trig),nWin*2+1);
    each.peak=nan(length(target),length(trig),nWin*2+1);
    for n=1:length(target)
        temp=coactTime.timestamp{target(n)};
        time=temp(temp>tRange(1) & temp<tRange(2))';
        val=coactTime.peakHeight{target(n)}(temp>tRange(1) & temp<tRange(2))';
        
        rateThis=nan(length(trig),nWin*2+1);
        for trigIdx=1:length(trig)
            cnt=histcounts(time,trig(trigIdx)+tBorder);
            rateThis(trigIdx,:)=cnt(2:end-1)/param.binSize;
            
            idx=cumsum(cnt);
            for m=1:length(cnt)-2
                strThis(trigIdx,m)=sum(val(idx(m)+1:idx(m+1)))/param.binSize;
                peakThis(trigIdx,m)=mean(val(idx(m)+1:idx(m+1)));
            end
        end
        rate(n,:)=mean(rateThis,1);
        strength(n,:)=mean(strThis,1);
        peak(n,:)=nanmean(peakThis,1);
        
        each.rate(n,:,:)=rateThis;
        each.strength(n,:,:)=strThis;
        each.peak(n,:,:)=peakThis;
    end

    
    %%
    
    evtTrigCoactChamber.rate=rate;
    evtTrigCoactChamber.peak=peak;
    evtTrigCoactChamber.strength=strength;
    evtTrigCoactChamber.each=each;
    evtTrigCoactChamber.time=tBin;
    evtTrigCoactChamber.tRange=tRange;    

    evtTrigCoactChamber.pairID=target;
    evtTrigCoactChamber.reacID=coactTime.reacID(target,:);
    evtTrigCoactChamber.region=coactTime.region(target,:);
    evtTrigCoactChamber.sigLevel=coactTime.sigLevel(target);
    evtTrigCoactChamber.tGap=coactTime.tGap(target);

    evtTrigCoactChamber.param=param;
    evtTrigCoactChamber.generator=mfilename;
    evtTrigCoactChamber.generatedate=datestr(now,'yyyy-mm-dd');

if ~strcmp(param.varName,'evtTrigCoactChamber')
    eval(sprintf('%s=evtTrigCoactChamber;',param.varName))
end
save(param.saveFile,param.varName,'-v7.3');









