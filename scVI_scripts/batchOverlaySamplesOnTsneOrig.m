topres='../results';
dirlocs={'simplenmf_ligernorm'}
allmerged='/mnt/dv/wid/projects5/Roy-singlecell/ke_work/engin_lab/monocle3_data_processing/monocle3_normalized_matrix.txt'
addpath('/mnt/dv/wid/projects5/Roy-singlecell/sr_work/multitask_matfact/matlab/');

perdatasetnormalize=1;
for d=1:length(dirlocs)
clear gnames cellnames datasets
[gnames,allcellnames,alldata]=loadDatasets_One(allmerged,perdatasetnormalize);
loadSampleSets
for k=10:5:20
	outdir=sprintf('%s/%s/k%d',topres,dirlocs{d},k);
	coord=load('tsne_coords_orig.mat')
	figfname=sprintf('%s/samples_on_tsne_orig.png',outdir);
	overlaySample(coord.mappedX_30,cellnames,sampleorder,figfname,2,5);

end
end
