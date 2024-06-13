topres='../scVI_application/results/';
dirlocs={'plots'};
scvifolder='/mnt/dv/wid/projects5/Roy-singlecell/ke_work/engin_lab/resequenced_data/scVI_application/results'
addpath('/mnt/dv/wid/projects5/Roy-singlecell/sr_work/multitask_matfact/matlab/');
addpath('/mnt/dv/wid/projects5/Roy-singlecell/sr_work/Flt-sne/FIt-SNE/');
%allmerged='/mnt/dv/wid/projects7/Roy-singlecell2/sr_work/scvi_tools/engin_results/scvi_engin_normalized_batch.txt'
%allmerged='/mnt/dv/wid/projects7/Roy-singlecell2/sr_work/scvi_tools/engin_results_batch/scvi_engin_normalized_batch.txt'
allmerged=sprintf('%s/resequenced_scVI_normalized_corrected_scaled_batch_t.txt',scvifolder);
alldata_corrected_flat_reorder=importdata(allmerged);
%[gnames,allcellnames,alldata]=loadDatasets_One(allmerged,perdatasetnormalize);
loadSampleSets
markers={'Ins1';'Ins2';'Gcg';'Sst';'Cd19';'Adgre1';'Il2ra';'Cd8a'};
%gnames=importdata('/mnt/dv/wid/projects7/Roy-singlecell2/sr_work/engin_project/data/monocle_corrected/genenames.txt');
gnames=importdata('/mnt/dv/wid/projects5/Roy-singlecell/ke_work/engin_lab/resequenced_data/scVI_application/input_data/genenames.txt');

for d=1:length(dirlocs)
for k=10:5:10%20
        outdir=sprintf('%s/%s/',topres,dirlocs{d});
        mkdir(outdir)
        coord=importdata(sprintf('%s/resequenced_scVI_mde_coord.txt',scvifolder));
        figfname=sprintf('%s/samples_on_scvi_sep.png',outdir);

        %Overlay samples on tsne
        overlaySample(coord,cellnames,sampleorder,figfname,2,5);
        continue
        %Overlay cluster on tsne
        fid=fopen(sprintf('../scVI_application/results/scvi_clusters.txt',outdir));
        cid=textscan(fid,'%s%d');
        fclose(fid);
        figfname=sprintf('%s/scvi_clusters_on_scvi.png',outdir);
        overlayClusterID_Tsne(coord,cid{2},6,5,figfname);
        figfname=sprintf('%s/scvi_clusters_on_scvi_oneplot.pdf',outdir);
        overlayClusterID_Tsne_oneplot(coord,cid{2},figfname);

        %fid=fopen(sprintf('%s/../cellclust_kmeans_norm_%d.txt',outdir,k));
        %cid=textscan(fid,'%s%d');
        %figfname=sprintf('%s/kmeans_on_tsne.png',outdir);
        %overlayClusterID_Tsne(coord,cid{2},k/5,5,figfname);

        figfname=sprintf('%s/markers_on_scvi_scaled.png',outdir);
        %plotMarkers_v2(markers,gnames,coord,alldata',2,4,figfname);
        plotMarkers_Engin(markers,gnames,coord,alldata_corrected_flat_reorder,2,4,figfname);
        %figfname=sprintf('%s/markers_on_tsne_bin.png',outdir);
        %plotMarkers_v2(markers,gnames,coord,alldata_bin',2,4,figfname);

end
end
