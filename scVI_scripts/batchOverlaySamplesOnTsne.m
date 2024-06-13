topres='../results';
dirlocs={'simplenmf_union_Tcells_ligernorm'}
dirlocs={'simplenmf_union_Tcells'}
dirlocs={'simplenmf_union_Tcells_allcells_ligernorm'}
dirlocs={'simplenmf_union_exp3_7_untransduced'};
dirlocs={'simplenmf_union_exp7_crx'};
%dirlocs={'simplenmf_union_exp7'}
dirlocs={'simplenmf_cowan_saha'};
dirlocs={'simplenmf_cowan_saha_organoid'};
dirlocs={'simplenmf_union_exp6_NRL1a_noundef'};
dirlocs={'simplenmf_union_exp6_NRL_Safe1';'simplenmf_union_exp6_NRL_Safe2'};
dirlocs={'simplenmf_union_exp6_NRL_Safe1_posonly';'simplenmf_union_exp6_NRL_Safe2'};
dirlocs={'simplenmf_union_exp6_UTF1';'simplenmf_union_exp6_UTF2'};
%dirlocs={'simplenmf_union_exp6_OrgUTF1'};

addpath('/mnt/dv/wid/projects5/Roy-singlecell/sr_work/multitask_matfact/matlab/');
addpath('/mnt/dv/wid/projects5/Roy-singlecell/sr_work/Flt-sne/FIt-SNE/');
%datadir='/mnt/dv/wid/projects5/Roy-singlecell/sahalab/data_U01/experiment7/depth_normalized/cas9_filtered_SR_union'
%samples={'count_Exp7_CRX1_neg';'count_Exp7_CRX1_pos';'count_Exp7_CRX2_neg';'count_Exp7_CRX2_pos';'count_Exp7_CRXSAFE_neg';'count_Exp7_CRXSAFE_pos';'count_Exp7_NRL1_neg';'count_Exp7_NRL1_pos';'count_Exp7_NRL2_neg';'count_Exp7_NRL2_pos';'count_Exp7_NRLSAFE_neg';'count_Exp7_NRLSAFE_pos';'count_Exp7_UTFORG';'count_Exp7_UTF'};
datadir='/mnt/dv/wid/projects5/Roy-singlecell/sahalab/data_U01/experiment3_7_cas9_filtered_SR_union'
%samples={'Untransduced';'count_Exp7_UTFORG';'count_Exp7_UTF'};
%datadir='/mnt/dv/wid/projects5/Roy-singlecell/sahalab/data_U01/experiment7/depth_normalized/cas9_filtered_SR_union'
%samples={'count_Exp7_CRX1_pos';'count_Exp7_CRX2_pos';'count_Exp7_CRX1_neg';'count_Exp7_CRX2_neg'};
%datadir='/mnt/dv/wid/projects5/Roy-singlecell/sahalab/integrate_cowan_saha/expression_matrices'
%samples={'Untransduced';'adult_human_retina_normal_fovea';'adult_human_retina_normal_periphery';'count_Exp7_UTF';'count_Exp7_UTFORG';'developed_human_retinal_organoid'};
datadir='/mnt/dv/wid/projects5/Roy-singlecell/sahalab/integrate_cowan_saha/developing_organoid/expression_matrices'
datadir='/mnt/dv/wid/projects5/Roy-singlecell/sahalab/data_U01/experiment6/depth_normalized'
samples={'Untransduced';'count_Exp7_UTF';'count_Exp7_UTFORG';'week_6';'week_12';'week_18';'week_24';'week_30';'week_38';'week_46'}
samples={'UTF1';'UTF2'};
samples={'cas9_pos_strict';'cas9_pos';'cas9_neg';'cas9_undef';};
%samples={'cas9_pos_strict';'cas9_pos';'cas9_neg';} %'cas9_undef';};
samplekey={'NRL_Safe1';'NRL_Safe2'};
samplekey={'UTF1';'UTF2'};
%samplekey={'OrgUTF1'};
samples={'cas9_pos_strict';'cas9_pos';} %'cas9_neg';'cas9_undef';};
samples={'cas9_neg';'cas9_undef';};
%[gnames,cellnames,datasets]=loadDatasets_func_v2(datadir,samples,perdatasetnormalize,'union');
perdatasetnormalize=1;
for d=1:length(dirlocs)
clear gnames cellnames datasets
[gnames,cellnames,datasets]=loadDatasets_func_Onesample(datadir,samples,perdatasetnormalize,samplekey{d});
for k=15:5:15
	outdir=sprintf('%s/%s/k%d',topres,dirlocs{d},k);
	coord=read_data(sprintf('%s/fast_tsne_nmf.dat',outdir));
	figfname=sprintf('%s/samples_on_tsne_sep.png',outdir);
	overlaySample(coord,cellnames,samples,figfname,2,5);

end
end
