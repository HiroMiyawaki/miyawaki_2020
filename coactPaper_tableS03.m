function coactPaper_tableS03
close all
height=15.5;
fh=initFig('height',height);

panel_01()
drawnow()

print(fh,'tableS03.pdf','-dpdf','-painters','-r300')

end

%%
function panel_01()
ica=poolVar('icaReacCCG_sig.mat');


ratList=fieldnames(ica);
%%
sesIdx=2;
hcIdx=3;

clear reg sig
reg={};
sig=[];
for ratIdx=1:length(ratList)
    rat=ratList{ratIdx};
    tempSig=ica.(rat)(sesIdx).nrem.significance(:,hcIdx);
    tempReg=ica.(rat)(sesIdx).region(ica.(rat)(sesIdx).pairID);
        
    idx=find(~strcmp(tempReg(:,1),tempReg(:,2)));
    reg=[reg;tempReg(idx,:)];
    sig=[sig;tempSig(idx)];
end
pairList={};
[pairList,~,pairIdx]=uniqueCellRows(reg);

sigList={'BLA','PrL L5'
         'vCA1','PrL L5'
         'vCA1','BLA'};
listOrder=[];
doFlip=[];
for n=1:size(sigList,1)
    idx=find(strcmp(pairList(:,1),sigList(n,1))&strcmp(pairList(:,2),sigList(n,2)));
    listOrder(end+1)=idx;
    doFlip(end+1)=false;
end

keyReg={'BLA','vCA1','PrL L5','vCA3','vSub','LA','CeA','PrL L5','PrL L2/3'};

for n=1:length(keyReg)
    idx=find(strcmp(pairList(:,1),keyReg{n}))';
    idx(ismember(idx,listOrder))=[];
    
    listOrder=[listOrder,idx];
    doFlip=[doFlip,false(size(idx))];
        
    idx=find(strcmp(pairList(:,2),keyReg{n}))';
    idx(ismember(idx,listOrder))=[];
    listOrder=[listOrder,idx];
    doFlip=[doFlip,true(size(idx))];
end

idx=find(~ismember(1:length(pairList),listOrder));
listOrder=[listOrder,idx];
doFlip=[doFlip,false(size(idx))];

%%
clear nEnsemble
nEnsemble.total=[];
nEnsemble.pos=[];
nEnsemble.neg=[];
nEnsemble.reg={};

for n=1:length(listOrder)
    if doFlip(n)==1
        nEnsemble.reg(n,:)=fliplr(pairList(listOrder(n),:));
    else
        nEnsemble.reg(n,:)=pairList(listOrder(n),:);
    end
    nEnsemble.total(n)=sum(pairIdx==listOrder(n));
    nEnsemble.pos(n)=sum(pairIdx==listOrder(n)&sig==1);
    nEnsemble.neg(n)=sum(pairIdx==listOrder(n)&sig==-1);
end

%%
clf


lineHeigth=5.5;
yMargin=0;
xMargin=2.5;
lineGap=0.3;

cellWidth=26;

height=lineHeigth*(length(nEnsemble.total)+2+lineGap)+yMargin*2;
width=cellWidth*4.5;

colNames={'Region pair','Number of pairs','Coupled pairs','Inverse-coupled pairs'};

subplotInMM(1,2,width,height)
hold on
ylim([0,height]);
xlim([0,width]);
set(gca,'YDir','reverse')

nLine=0;


plot([0,width],[0,0]+lineHeigth*(nLine+lineGap)+yMargin,'k-','LineWidth',0.75)
nLine=nLine+1;
for n=1:length(colNames)
    text((cellWidth+xMargin)*(n-1)+xMargin,lineHeigth*nLine+yMargin,colNames{n},...
        'fontsize',7,'horizontalALign','left','verticalAlign','baseline')
end
plot([0,width],lineHeigth*(nLine+lineGap)+yMargin+[0,0],'k-','LineWidth',0.75)

for m=1:length(nEnsemble.total)
nLine=nLine+1;
    for n=1:length(colNames)
        txt=[];
        switch n
            case 1
                txt=join(nEnsemble.reg(m,:), ' - ');
            case 2
                txt=num2str(nEnsemble.total(m))
            case 3
                txt=sprintf('%d (%0.1f%%)',nEnsemble.pos(m),nEnsemble.pos(m)/nEnsemble.total(m)*100)
            case 4
                txt=sprintf('%d (%0.1f%%)',nEnsemble.neg(m),nEnsemble.neg(m)/nEnsemble.total(m)*100)
        end
        text((cellWidth+xMargin)*(n-1)+xMargin,lineHeigth*nLine+yMargin,txt,...
            'fontsize',7,'horizontalALign','left','verticalAlign','baseline')
    end
end

plot([0,width],lineHeigth*(nLine+lineGap)+yMargin+[0,0],'k-','LineWidth',0.75)
axis off
end
