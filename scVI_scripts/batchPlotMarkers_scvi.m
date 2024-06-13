topres='../scVI_application/';
dirlocs={'results'}
markers={'Ins1';'Ins2';'Gcg';'Sst';'Cd19';'Adgre1';'Il2ra';'Cd8a'}
%markers={'Pou5f1';'Esrrb';'Hes1';'Gmeb1';'Nr5a2';'BFP'};
addpath('/mnt/dv/wid/projects5/Roy-singlecell/sr_work/multitask_matfact/matlab/');
addpath('/mnt/dv/wid/projects5/Roy-singlecell/sr_work/Flt-sne/FIt-SNE/');
alldata_bin=alldata;
nz=find(alldata_bin>0);
alldata_bin(nz)=1;
for d=1:1 %length(dirlocs)
%for k=10:10:30
%for k=10:10:30
for k=10
	outdir=sprintf('%s/%s',topres,dirlocs{d});
	coord=importdata(sprintf('%s/resequenced_scVI_mde_coord.txt',outdir));
	figfname=sprintf('%s/markers_on_tsne2.pdf',outdir);
	plotMarkers_Engin_v3(markers,gnames,coord,alldata',3,5,figfname);
	%figfname=sprintf('%s/markers_on_tsne_bin.png',outdir);
	%plotMarkers_v2(markers,gnames,coord,alldata_bin',3,5,figfname);
end
end
