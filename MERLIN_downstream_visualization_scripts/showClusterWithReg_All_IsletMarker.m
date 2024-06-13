function csize=showClusterWithReg_Islet_All(clust,clustreg,alldata1,allregdata,minsize,samplelocs,samplenames,figfname,isletexp,isletmarkers)
gclust=clust.data;
gnames=clust.textdata;
regclust=clustreg.data;
regnames=clustreg.textdata;

cids=unique(gclust);
mattoshow=[];
csize=histc(gclust,cids);
toshowcids=find(csize>minsize);
fprintf('Found %d modules of size at least %d\n',length(toshowcids),minsize);
%colormap(redbluecmap);
v1=(0:0.01:1)';
v2=(1:-0.01:0)';
v3=ones(101,1);
c=[v1,v1,v3;v3,v2,v2];
colormap(c);
clusterparts=[];
samplepos=[];
for i=1:length(toshowcids)
moduleofinterest=cids(toshowcids(i));
ids=find(gclust==moduleofinterest);
if(length(ids)<minsize)
	fprintf('Skiping module:: %d size=%d\n',moduleofinterest,length(ids));
	continue;
end
hold off
subplot(1,1,1);
%%module expression
cmat=alldata1(ids,:);
cmat=cmat-repmat(mean(cmat,2),1,size(cmat,2));
%%module regulators expression
regids=find(regclust==moduleofinterest);
regmat=allregdata(regids,:);
regmat=regmat-repmat(mean(regmat,2),1,size(regmat,2));
toshowmat=[cmat;ones(1,size(regmat,2));regmat];
%%Islet marker expression
isletexp=isletexp-repmat(mean(isletexp,2),1,size(cmat,2));
toshowmat=[toshowmat;ones(1,size(regmat,2));isletexp];
q=2;

imagesc(toshowmat,[-0.5 0.5]);
set(gca,'yticklabels',[gnames(ids);' ';regnames(regids);' ';isletmarkers],'fontsize',6);
yticks([1:size(toshowmat,1)]);
colorbar;
cellcnt=0;
for c=1:length(samplelocs)
        cnt=samplelocs(c);
        cellcnt=cellcnt+cnt;
        line([cellcnt cellcnt],[0 size(toshowmat,1)],'color',[0 0 0],'linewidth',2);
	samplepos(c)=cellcnt;
end
set(gca,'xticklabels',strrep(samplenames,'_','-'),'fontsize',6);
xticks(samplepos);
xtickangle(45);
mattoshow=[mattoshow;cmat]; 
clusterparts=[clusterparts length(ids)];
height=size(toshowmat,1)*0.1;
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 height+1],'PaperSize',[6 height+1]);
saveas(gcf,sprintf('%s/C%d.png',figfname,cids(toshowcids(i))),'png');
end

subplot(1,1,1);
[ig,order]=sort(gclust);
imagesc(mattoshow,[-0.5 0.5]);
cellcnt=0;
for c=1:length(samplelocs)
        cnt=samplelocs(c);
        cellcnt=cellcnt+cnt;
        line([cellcnt cellcnt],[0 size(alldata1,1)],'color',[0 0 0],'linewidth',2);
end
gcnt=0;
gcntset=[];
for c=1:length(clusterparts)
	cnt=clusterparts(c);
	gcnt=gcnt+cnt;
        line([0 size(alldata1,2)],[gcnt gcnt],'color',[0 0 0],'linewidth',0.1);
	gcntset=[gcntset gcnt];
end
%set(gca,'xticklabels',strrep(samplenames,'_','-'),'xticks',samplepos);
set(gca,'xticklabels',strrep(samplenames,'_','-'),'fontsize',6);
xticks(samplepos);
xtickangle(45);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 6 10],'PaperSize',[6 10]);
set(gca,'yticklabels',cids(toshowcids));
yticks(gcntset);
colorbar;
saveas(gcf,sprintf('%s/allcids_5-9_temp.png',figfname),'png');
