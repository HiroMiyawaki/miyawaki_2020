function fear_testSignificance_reacCCG(basename,varargin)
% basename='~/data/Fear/triple/achel180320/achel180320';
load([basename '.basicMetaData.mat'])
fprintf('%s start %s with data of %s\n',datestr(now),mfilename,basicMetaData.SessionName)

param.maxDelta=100; %ms
param.redo=false; %ms
param=parseParameters(param,varargin);
%%

varList={};
if exist([basicMetaData.AnalysesName '-instReacCCG_sh.mat'],'file')
    varList{end+1}='instReacCCG_sh';
end
if exist([basicMetaData.AnalysesName '-icaReacCCG_sh.mat'],'file')
    varList{end+1}='icaReacCCG_sh';
end

if exist([basicMetaData.AnalysesName '-instReacCCG_exSWR_sh.mat'],'file')
    varList{end+1}='instReacCCG_exSWR_sh';
end

if exist([basicMetaData.AnalysesName '-icaReacCCG_exSWR_sh.mat'],'file')
    varList{end+1}='icaReacCCG_exSWR_sh';
end

if exist([basicMetaData.AnalysesName '-instReacZNCCG.mat'],'file')
    varList{end+1}='instReacZNCCG';
end
if exist([basicMetaData.AnalysesName '-icaReacZNCCG.mat'],'file')
    varList{end+1}='icaReacZNCCG';
end

if exist([basicMetaData.AnalysesName '-instReacZNCCG_exSWR.mat'],'file')
    varList{end+1}='instReacZNCCG_exSWR';
end

if exist([basicMetaData.AnalysesName '-icaReacZNCCG_exSWR.mat'],'file')
    varList{end+1}='icaReacZNCCG_exSWR';
end

if exist([basicMetaData.AnalysesName '-instReacCCGchamber_sh.mat'],'file')
    varList{end+1}='instReacCCGchamber_sh';
end
if exist([basicMetaData.AnalysesName '-instReacZNCCGchamber.mat'],'file')
    varList{end+1}='instReacZNCCGchamber';
end

if exist([basicMetaData.AnalysesName '-icaReacCCGchamber_sh.mat'],'file')
    varList{end+1}='icaReacCCGchamber_sh';
end

if exist([basicMetaData.AnalysesName '-icaReacZNCCGchamber.mat'],'file')
    varList{end+1}='icaReacZNCCGchamber';
end


if exist([basicMetaData.AnalysesName '-instReacCCG_exKcomplex_sh.mat'],'file')
    varList{end+1}='instReacCCG_exKcomplex_sh';
end

if exist([basicMetaData.AnalysesName '-icaReacCCG_exKcomplex_sh.mat'],'file')
    varList{end+1}='icaReacCCG_exKcomplex_sh';
end


if exist([basicMetaData.AnalysesName '-icaReacZNCCG_exKcomplex.mat'],'file')
    varList{end+1}='icaReacZNCCG_exKcomplex';
end

if exist([basicMetaData.AnalysesName '-instReacZNCCG_exKcomplex.mat'],'file')
    varList{end+1}='instReacZNCCG_exKcomplex';
end

if exist([basicMetaData.AnalysesName '-instReacZNCCGchamberCue.mat'],'file')
     varList{end+1}='instReacZNCCGchamberCue';
end
    
if exist([basicMetaData.AnalysesName '-icaReacZNCCGchamberCue.mat'],'file')
     varList{end+1}='icaReacZNCCGchamberCue';
end

if exist([basicMetaData.AnalysesName '-icaReacZNCCG_exHFObaseCond.mat'],'file')
    varList{end+1}='icaReacZNCCG_exHFObaseCond';
end

%%
if isempty(varList)
    warning('Nothing to be processed: exit')
    return
end

%%
maxDelta=param.maxDelta;

for varIdx=1:length(varList);
    
    if strcmp(varList{varIdx}(end-2:end),'_sh')
        resName=[varList{varIdx}(1:end-2) 'sig'];
    else
        resName=[varList{varIdx} '_sig'];
    end
    
    if  ~param.redo && exist([basicMetaData.AnalysesName '-' resName '.mat'],'file') 
        fprintf('%s already exist!\n',resName)
        continue
    end
    load([basicMetaData.AnalysesName '-' varList{varIdx} '.mat'])
    target=eval(varList{varIdx});

    
    if ismember(varList{varIdx},{'instReacCCG_sh','icaReacCCG_sh','instReacZNCCG','icaReacZNCCG'})
        behList={'entire','wake','nrem','rem'};
    elseif ismember(varList{varIdx},{'instReacCCG_exSWR_sh','icaReacCCG_exSWR_sh','instReacZNCCG_exSWR','icaReacZNCCG_exSWR'})
        behList={'exRipple','exSpindle','exBoth'};
    elseif ismember(varList{varIdx},{'instReacCCG_exKcomplex_sh','icaReacCCG_exKcomplex_sh','instReacZNCCG_exKcomplex','icaReacZNCCG_exKcomplex'})
        behList={'exKcomp','exKgap'};
    elseif ismember(varList{varIdx},{'icaReacZNCCG_exHFObaseCond'})
        behList={'exHFO'};
    elseif ismember(varList{varIdx},{'instReacCCGchamber_sh','instReacZNCCGchamber','icaReacCCGchamber_sh','icaReacZNCCGchamber','instReacZNCCGchamberCue','icaReacZNCCGchamberCue'})
        behList={};
    else
        warning('Unknown var name: %s',varList{varIdx})
        continue
    end
    fprintf('    %s processing %s\n',datestr(now),varList{varIdx})
    
    clear res
    for tempIdx=1:length(target);
        tBinSize=target(tempIdx).tBinSize*1000;
        
        if ~isempty(behList)
            centerBin=(size(target(tempIdx).(behList{1}).real.ccg,2)+1)/2;
            peakRange=ceil(maxDelta/tBinSize);
            peakRange=(-peakRange:peakRange)+centerBin;
            for behIdx=1:length(behList);
                beh=behList{behIdx};
                delta=target(tempIdx).(beh).real.ccg-target(tempIdx).(beh).shuffle.mean;
                
                [~,pos]=max(abs(delta(:,peakRange,:)),[],2);
                pos=pos+(peakRange(1)-1);
                
                
                for pIdx=1:size(target(tempIdx).(beh).real.ccg,1)
                    for tarIdx=1:size(target(tempIdx).(beh).real.ccg,3)
                        val=target(tempIdx).(beh).real.ccg(pIdx,pos(pIdx,1,tarIdx),tarIdx);
                        res(tempIdx).(beh).peakTime(pIdx,tarIdx)=(pos(pIdx,1,tarIdx)-centerBin)*tBinSize;
                        res(tempIdx).(beh).peakValue(pIdx,tarIdx)=delta(pIdx,pos(pIdx,1,tarIdx),tarIdx);
                        res(tempIdx).(beh).significance(pIdx,tarIdx)=(val>target(tempIdx).(beh).shuffle.global99(pIdx,2,tarIdx)) - (val<target(tempIdx).(beh).shuffle.global99(pIdx,1,tarIdx));
                        res(tempIdx).(beh).significance5(pIdx,tarIdx)=(val>target(tempIdx).(beh).shuffle.global95(pIdx,2,tarIdx)) - (val<target(tempIdx).(beh).shuffle.global95(pIdx,1,tarIdx));
                    end
                end
            end
        else
            centerBin=(size(target(tempIdx).real.ccg,2)+1)/2;
            peakRange=ceil(maxDelta/tBinSize);
            peakRange=(-peakRange:peakRange)+centerBin;
            delta=target(tempIdx).real.ccg-target(tempIdx).shuffle.mean;
            
            [~,pos]=max(abs(delta(:,peakRange,:)),[],2);
            pos=pos+(peakRange(1)-1);
            
            
            for pIdx=1:size(target(tempIdx).real.ccg,1)
                for tarIdx=1:size(target(tempIdx).real.ccg,3)
                    val=target(tempIdx).real.ccg(pIdx,pos(pIdx,1,tarIdx),tarIdx);
                    res(tempIdx).peakTime(pIdx,tarIdx)=(pos(pIdx,1,tarIdx)-centerBin)*tBinSize;
                    res(tempIdx).peakValue(pIdx,tarIdx)=delta(pIdx,pos(pIdx,1,tarIdx),tarIdx);
                    res(tempIdx).significance(pIdx,tarIdx)=(val>target(tempIdx).shuffle.global99(pIdx,2,tarIdx)) - (val<target(tempIdx).shuffle.global99(pIdx,1,tarIdx));
                    res(tempIdx).significance5(pIdx,tarIdx)=(val>target(tempIdx).shuffle.global95(pIdx,2,tarIdx)) - (val<target(tempIdx).shuffle.global95(pIdx,1,tarIdx));
                end
            end
        end
        res(tempIdx).pairID=target(tempIdx).pairID;
        res(tempIdx).region=target(tempIdx).region;
        res(tempIdx).instReacID=target(tempIdx).instReacID;
        res(tempIdx).template=target(tempIdx).template;
        res(tempIdx).generator=mfilename;
        res(tempIdx).generatedate=datestr(now,'yyyy-mm-dd');
        res(tempIdx).param=param;
        
    end
    eval(sprintf('%s=res;',resName))
    
        save([basicMetaData.AnalysesName '-' resName '.mat'],resName,'-v7.3')
end



