function fear_homecage_spkCCG_1ms(basename)
load([basename '.basicMetadata.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.AnalysesName '-spkCCG.mat'])
load([basicMetaData.AnalysesName '-okUnit.cellinfo.mat'])

%%
[reg,regList]=relabel_region(okUnitInfo.region,'minCellNum',5);
regList(strcmpi(regList,'other'))=[];

%%
regPair=[];
for n=1:length(regList)
    for m=1:n
        regPair(end+1,:)=[m,n];
    end
end

[~,order]=sortrows([diff(regPair,1,2)~=0,regPair(:,1)])
regPair=regPair(order,:)


for pIdx=1:size(regPair,1)
    preRegName=regList{regPair(pIdx,1)};
    postRegName=regList{regPair(pIdx,2)};

    preReg=find(strcmp(reg,preRegName));
    postReg=find(strcmp(reg,postRegName));



    typeList=fieldnames(spkCCG.homecage(1));
    for typeIdx=1:length(typeList)
        type=typeList{typeIdx};

        nBin=round(spkCCG.param.halfWindow(1)/spkCCG.param.binSize(1))*2+1;
        subset.(type)=zeros(nBin,length(preReg),length(postReg));
        for n=1:5
            if any(spkCCG.homecage(n).(type)(1).nSpk>0)
                subset.(type)=subset.(type)+spkCCG.homecage(n).(type)(1).smoothed(:,preReg,postReg);
            end
        end
        t=spkCCG.homecage(1).entire(1).t*1000;

        regCCG(pIdx).(type)=[];

        for n=1:size(preReg,2)
            for m=1:size(postReg,2)
                if preReg(n)==postReg(m)
                else
                    regCCG(pIdx).(type)(end+1,:)=subset.(type)(:,n,m);
                end
            end
        end
    end
    
    
    
end
%%
close all
fh=initFig4A4('fontsize',7)
tRange=20;
zRange=50;
for pIdx=1:size(regPair,1)
    preRegName=regList{regPair(pIdx,1)};
    postRegName=regList{regPair(pIdx,2)};
    
    temp=zscore(regCCG(pIdx).entire')';
    [~,order]=sort(sum(temp(:,(size(regCCG(pIdx).entire,2)+1)/2+(-1:1)),2),'descend');
    zT=abs(t)<=zRange;
    showT=t(zT);
    showTrange=abs(showT)<=tRange;

    for typeIdx=1:length(typeList)
        subplot2(size(regPair,1),6,pIdx,typeIdx)
        type=typeList{typeIdx};
        temp=zscore(regCCG(pIdx).(type)(:,zT)')';
        imagesc(showT(showTrange),[],temp(order,showTrange))
        title(type)
        box off
        ax=fixAxis;
        if typeIdx==1
            ylabel(sprintf('%s-%s\ncell pair',preRegName,postRegName))
        end
        if pIdx==size(regPair,1)
            xlabel('time (ms)')
        end
    end
end
addScriptName(mfilename,true)
textInMM(10,15,basicMetaData.SessionName,'fontsize',12)
%%
print(fh,[basicMetaData.AnalysesPdf ,'-fear_homecage_spkCCG_1ms.pdf'],'-dpdf')
