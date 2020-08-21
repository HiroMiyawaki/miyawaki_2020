function coactPaper_tableS01
close all
height=10;
width=2.54*6.5;
fh=initFig('height',height);


set(fh,'position',[0,20,(width-0),(height-0)])
set(fh, 'Units','centimeters')
set(fh,'PaperSize',[width,height])

panel_01()
drawnow()
print(fh,'tableS01.pdf','-dpdf','-painters','-r300')

end
%%
function panel_01()
okUnit=poolVar('okUnit.cellinfo.mat');
ngUnit=poolVar('ngUnit.cellinfo.mat');
instReac=poolVar('instReacInfo.mat');


ratList=fieldnames(okUnit);
%%
regList={'vCA1','vCA3','vSub','BLA','LA','CeA','PrL L2/3','PrL L5','Pir','STIA','Other'};
for ratIdx=1:length(ratList)
    temp=relabel_region(okUnit.(ratList{ratIdx}).region,'minCellNum',0);
    temp(~ismember(temp,regList))={'Other'};
    okUnit.(ratList{ratIdx}).region=temp;
    
    temp=relabel_region(ngUnit.(ratList{ratIdx}).region,'minCellNum',0);
    temp(~ismember(temp,regList))={'Other'};
    ngUnit.(ratList{ratIdx}).region=temp;
end

%%

%%
%%
clear nCell
nCell.ok=[];
nCell.ng=[];
nCell.ex=[];
nCell.inh=[];
nCell.nc=[];
for regIdx=1:length(regList)
    for ratIdx=1:length(ratList)
        nCell.ok(ratIdx,regIdx)=sum(strcmpi(okUnit.(ratList{ratIdx}).region,regList{regIdx}));
        nCell.ex(ratIdx,regIdx)=sum(strcmpi(okUnit.(ratList{ratIdx}).region,regList{regIdx})&okUnit.(ratList{ratIdx}).cellType.type==1);
        nCell.inh(ratIdx,regIdx)=sum(strcmpi(okUnit.(ratList{ratIdx}).region,regList{regIdx})&okUnit.(ratList{ratIdx}).cellType.type==-1);
        nCell.nc(ratIdx,regIdx)=sum(strcmpi(okUnit.(ratList{ratIdx}).region,regList{regIdx})&okUnit.(ratList{ratIdx}).cellType.type==0);
        nCell.ng(ratIdx,regIdx)=sum(strcmpi(ngUnit.(ratList{ratIdx}).region,regList{regIdx}));
    end
end
%%

sesIdx=2;
temp={};
phi=[];
rat=[];
for ratIdx=1:length(ratList)
    ratName=ratList{ratIdx};
    phi=instReac.(ratName)(sesIdx).phi;
    temp=relabel_region(instReac.(ratName)(sesIdx).region(phi>1),'minCellNum',0);
    temp=relabel_region(temp,'minCellNum',0);
    for regIdx=1:length(regList)
        nCell.ica(ratIdx,regIdx)=sum(strcmpi(temp,regList{regIdx}));
    end
end
%%
mesList=fieldnames(nCell)
for n=1:length(mesList)
    nCell.(mesList{n})(:,end+1)=sum(nCell.(mesList{n}),2);
    nCell.(mesList{n})(end+1,:)=sum(nCell.(mesList{n}),1);
end

nCell.total=nCell.ng+nCell.ok;
regList(end+1)={'Total'};
ratList{end+1}='Total';
%%
headCap=@(x) [upper(x(1)) x(2:end)];
%%
%%
width=163;

ratNameWidth=15;
lineHeigth=5;
yMargin=0.5;
xMargin=1;
lineGap=0.6;

cellWidth=(width-ratNameWidth-xMargin)/(length(regList)-1)-xMargin;

height=lineHeigth*(2*(length(ratList)+3)+2)+yMargin*2;

subplotInMM(1,2,width,height)
hold on
ylim([0,height]);
xlim([0,width]);
set(gca,'YDir','reverse')

nLine=0;


plot([0,width],[0,0]+yMargin,'k-','LineWidth',0.75)
for type=1

nLine=nLine+1;
text(ratNameWidth+(cellWidth+xMargin)*(0+length(regList)-1)/2,lineHeigth*nLine,'Brain region',...
    'fontsize',7,'horizontalALign','center','verticalAlign','baseline')
plot(ratNameWidth+(cellWidth+xMargin)*[0,length(regList)-1],lineHeigth*(nLine+lineGap)-1.5+[0,0],'k-','LineWidth',0.75)

nLine=nLine+1;
yPos=lineHeigth*nLine;
for regIdx=0:length(regList)-1
    if regIdx==0
        xPos=xMargin;
        text(xPos,yPos,'Rat name','horizontalALign','left','fontsize',7,'verticalAlign','baseline')
    else
        xPos=xMargin+ratNameWidth+(regIdx-1)*(cellWidth+xMargin);
        text(xPos,yPos,regList{regIdx},'horizontalALign','left','fontsize',7,'verticalAlign','baseline')
    end
end
plot([0,width-xMargin],[0,0]+yMargin+(lineHeigth*nLine+lineGap),'k-','LineWidth',0.75)

for ratIdx=1:length(ratList)
    nLine=nLine+1;
    yPos=lineHeigth*nLine;
    xPos=xMargin;
    if ratIdx<length(ratList)
        text(xPos,yPos,headCap(ratList{ratIdx}(1:end-6)),'horizontalALign','left','fontsize',7)    
    else
        plot([xMargin,width-xMargin],lineHeigth*(nLine-1+lineGap)+[0,0],'k-','LineWidth',0.75)
        text(xPos,yPos,headCap(ratList{ratIdx}),'horizontalALign','left','fontsize',7)    
    end
    for regIdx=1:length(regList)-1
        xPos=xMargin+ratNameWidth+(regIdx-1)*(cellWidth+xMargin);
        if nCell.ok(ratIdx,regIdx)+nCell.ng(ratIdx,regIdx)>0
            if type==1
                text(xPos,yPos,sprintf('%d, %d, %d',...
                    nCell.ex(ratIdx,regIdx),nCell.inh(ratIdx,regIdx),nCell.nc(ratIdx,regIdx)),...
                    'horizontalALign','left','fontsize',7)
            else
                text(xPos,yPos,sprintf('%d',...
                    nCell.ica(ratIdx,regIdx)),...
                    'horizontalALign','left','fontsize',7)
            end
        else
            text(xPos+0.5,yPos,'--',...
                'horizontalALign','center','fontsize',7,'Interpreter','latex')
        end
    end
end
plot([0,width],lineHeigth*(nLine+lineGap)+[0,0],'k-','LineWidth',0.75)
nLine=nLine+1;
end
axis off
end

