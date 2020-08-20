function fear_amyGammaTrigIcaCoactStrength(basename,varargin)
%%
param.templateIdx=2; %id of template behavior session
param.targetIdx=3; %id of homecage session for significant selection
param.beh='nrem'; %name of behaviror type for siginificant selection : 'nrem' 'rem' 'wake' or 'entire'
param.tWin=5;

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)
%%
param=parseParameters(param,varargin);
%%
load([basicMetaData.AnalysesName '-icaReac.mat'])
load([basicMetaData.AnalysesName '-icaReacZNCCG_sig.mat'])

load([basicMetaData.Basename '.sessions.events.mat'])
load([basicMetaData.Basename '.cues.events.mat'])
load([basicMetaData.Basename '.SleepState.states.mat'])

load([basicMetaData.Basename '.amyLowGamma.events.mat'])
load([basicMetaData.Basename '.amyHighGamma.events.mat'])
%%
slp=relabel_ma2sleep(SleepState.MECE.timestamps);

%%
tBinSize=icaReac(1).param.tBinSize*1e3;
templateIdx=param.templateIdx;
targetIdx=param.targetIdx;
beh=param.beh;
isSig=icaReacZNCCG_sig(templateIdx).(beh).significance(:,targetIdx) + ...
    icaReacZNCCG_sig(templateIdx).(beh).significance5(:,targetIdx);
reg=icaReacZNCCG_sig(templateIdx).region;




%%
fprintf('%s getting inst coact strength \n',datestr(now))

across=arrayfun(@(x,y) ~strcmpi(reg{x},reg{y}), icaReacZNCCG_sig(templateIdx).pairID(:,1),icaReacZNCCG_sig(templateIdx).pairID(:,2));

target=find(across);

targetPair=icaReacZNCCG_sig(templateIdx).pairID(target,:);
reacID=icaReacZNCCG_sig(templateIdx).instReacID(targetPair);
regPair=reg(targetPair);

sigLevel=zeros(size(target));
sigLevel(isSig(target)==2)=1;
sigLevel(isSig(target)==1)=5;
sigLevel(isSig(target)==-2)=-1;
sigLevel(isSig(target)==-1)=-5;


gap=icaReacZNCCG_sig(templateIdx).(beh).peakTime(target,targetIdx)/tBinSize;

icaCocatStrength=zeros(length(target),size(icaReac(templateIdx).strength,2));
for idx=1:length(target)
    
    x=icaReac(templateIdx).strength(reacID(idx,1),:);
    y=icaReac(templateIdx).strength(reacID(idx,2),:);
    
    x=zscore(x);
    y=zscore(y);
    
    if gap(idx)<0
        y=[y(1-gap(idx):end),zeros(1,-gap(idx))];
    else
        y=[zeros(1,gap(idx)),y(1:end-gap(idx))];
    end
    icaCocatStrength(idx,:)=x.*y;
    
end


tBin=(1:size(icaCocatStrength,2))*tBinSize;

%%
tWin=param.tWin;
nWin=ceil(tWin*1e3/tBinSize);

for tRangeType=1:2
    clear res
    if tRangeType==1
        tRange=sessions.homecage(3,:);
        typeMax=2;
        behIdx=3; %nrem
        fileName='amyGammaTrigCoact-postNREM.mat';
        varName='amyGammaTrigCoact';
        fprintf('%s getting peri event triggered average in post-NREM \n',datestr(now))
    else
        tMin=sessions.timestamps(4,1);
        tMax=sessions.timestamps(4,2);
        
        tMin=min(cues.timestamps.Pip(cues.timestamps.Pip(:,1)>tMin,1));
        
        tRange=[tMin,tMax];
        typeMax=2;
        behIdx=1; %awake
        fileName='amyGammaTrigCoact-cueRet.mat';
        varName='amyGammaTrigCoact';
        fprintf('%s getting peri event triggered average in cue-ret \n',datestr(now))
    end
    
    beh=slp(slp(:,2)>tRange(1)&slp(:,1)<tRange(2),:);
    if beh(1,1)<tRange(1);beh(1,1)=tRange(1);end
    if beh(end,2)>tRange(1);beh(end,2)=tRange(2);end
    
    beh=beh(beh(:,3)==behIdx,1:2);
    
    for trigType=1:typeMax
        if trigType==1
                trigT=amyLowGamma.peaks.timestamps;
                gmParam=amyLowGamma.param;
        elseif trigType==2
                trigT=amyHighGamma.peaks.timestamps;
                gmParam=amyHighGamma.param;
        end
        
        trigT=trigT(trigT>tRange(1) & trigT<tRange(2));
        
        trigT=trigT(any(trigT>beh(:,1)'&trigT<beh(:,2)',2));
        
        trigF=round(trigT*1e3/tBinSize);
        
        res.peth(trigType).avg=nan(size(icaCocatStrength,1),nWin*2+1);
        if ~isempty(trigF)
            for idx=1:size(icaCocatStrength,1)
                temp=icaCocatStrength(idx,:);
                res.peth(trigType).avg(idx,:)=mean(temp((-nWin:nWin)+trigF),1);
            end
        end
        res.peth(trigType).nTrig=length(trigT);
        res.peth(trigType).param=gmParam;
    end
    res.pairID=target;
    res.reacID=reacID;
    res.region=regPair;
    res.tGap=gap;
    res.sigLevel=sigLevel;
    res.param=param;
    res.generator=mfilename;
    res.generatedate=datestr(now,'yyyy-mm-dd');
    
    eval(sprintf('%s=res;',varName));
    
    fprintf('%s savig results \n',datestr(now))
    save([basicMetaData.AnalysesName '-' fileName],varName,'-v7.3')
end












