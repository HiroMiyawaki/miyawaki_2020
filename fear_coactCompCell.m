function fear_coactCompCell(basename,varargin)

load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

load([basicMetaData.AnalysesName '-icaReacInfo.mat'])
load([basicMetaData.AnalysesName '-icaReacPartner.mat'])
load([basicMetaData.AnalysesName '-okUnit.cellinfo.mat'])

load([basicMetaData.AnalysesName '-instReacInfo.mat'])
load([basicMetaData.AnalysesName '-instReacPartner.mat'])
%%
param.threshold=0.3;

param=parseParameters(param,varargin);
%%

[cellReg,regList]=relabel_region(okUnitInfo.region,'minCellNum',0);
nCell=length(okUnitInfo.region);
nReg=length(regList);

for reacType=1:2
    if reacType==1
        reacPartner=icaReacPartner;
        reacInfo=icaReacInfo;
        reacName='ica';
    else
        reacPartner=instReacPartner;
        reacInfo=instReacInfo;
        reacName='pca';
    end
    
    for tempSes=1:length(reacPartner);
        for matchHC=1:length(reacPartner(tempSes).partner);
            behList=fieldnames(reacPartner(tempSes).partner(matchHC));
            for bIdx=1:length(behList)
                beh=behList{bIdx};
                partnerList=reacPartner(tempSes).partner(matchHC).(beh).pos;
                
                for idx=1:length(partnerList)
                    partnerReg=(reacPartner(tempSes).region(partnerList{idx}));
                    for rIdx=1:length(regList)
                        nPair(idx,rIdx)=sum(strcmp(partnerReg,regList{rIdx}));
                    end
                end
                
                conn=zeros(nCell,nReg);
                for idx=1:length(partnerList)
                    cmpIdx=reacPartner(tempSes).instReacID(idx);
                    cellID=find(strcmp(cellReg,reacInfo(tempSes).region{cmpIdx}));
                    
                    if length(cellID)~=length(reacInfo(tempSes).weigth{cmpIdx})
                        error('number of cells and weigths are not matched')
                    end
                    
                    conn(cellID(abs(reacInfo(tempSes).weigth{cmpIdx})>param.threshold),nPair(idx,:)>0)=...
                        conn(cellID(abs(reacInfo(tempSes).weigth{cmpIdx})>param.threshold),nPair(idx,:)>0)+1;
                end
                
                coactCompCell.(reacName)(tempSes).homecage(matchHC).(beh)=conn;
                coactCompCell.(reacName)(tempSes).reacID=reacPartner(tempSes).instReacID;
                coactCompCell.(reacName)(tempSes).template=reacInfo(tempSes).tempName;
                
            end
        end
    end
end
coactCompCell.region=regList;
coactCompCell.generator=mfilename;
coactCompCell.generatedate=datestr(now,'yyyy-mm-dd');

%%
save([basicMetaData.AnalysesName '-coactCompCell.mat'],'coactCompCell','-v7.3')







