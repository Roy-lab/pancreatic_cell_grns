cluster='beta';
sampleloc{1}=[446,300,796,548,718,1986,2775];
samplenames={'KO_ire1_6680'; 'KO_ire1_6681'; 'WT_ire1_6683';'KO_xbp1_2117';'KO_xbp1_2119';'WT_xbp1_2120';'WT_xbp1_2118'};


DIRNAME=sprintf('/mnt/dv/wid/projects5/Roy-singlecell/ke_work/engin_lab/resequenced_data/MERLIN_on_scVI_corrected_reseq_data');
CLUSTERDIR=sprintf('%s/run_merlin_beta/modules',DIRNAME);
CLUSTERDIR_REG=sprintf('%s/run_merlin_beta/list_of_terms',DIRNAME);
SRCLUSTERDIR=sprintf('/mnt/dv/wid/projects5/Roy-singlecell/ke_work/engin_lab/resequenced_data/MERLIN_on_scVI_corrected_reseq_data/result');
EXPDIR='/mnt/dv/wid/projects5/Roy-singlecell/ke_work/engin_lab/resequenced_data/MERLIN_on_scVI_corrected_reseq_data/beta_cells_C0_3_data'
EXPFILE=sprintf('%s/beta_cells_C0_3_scVI_corrected_scaled_sqrt_tr.txt',EXPDIR);
EXPFILE
alldata=importdata(EXPFILE);
gnames=alldata.textdata(2:end,1);
alldata=alldata.data;
figure;
t=0.4   
        CASSIGN=sprintf('%s/module_%s.%.1f_geneset.txt',CLUSTERDIR,cluster,t);
        cid=importdata(CASSIGN);
        markers=cid.textdata;
        [gid,actualid]=getGeneIDs(gnames,markers);
        expdata=alldata(gid,:);
        fprintf('Found %d genes\n',length(gid));

        pseudobulkexp=getPseudoBulkExpPerSample(expdata,sampleloc{1},samplenames);

        wtsamples=[3,6,7];
        xbpsamples=[4,5];
        ire1samples=[1,2];
        fprintf('Now getting xbp1 modules\n')
        [xbp_modules_up, xbp_modules_down]=getDiffModules(cid,pseudobulkexp,xbpsamples,wtsamples);
        fprintf('Now getting xbp1 modules vs all other \n')
        [xbp_modules_up2, xbp_modules_down2]=getDiffModules(cid,pseudobulkexp,xbpsamples,[wtsamples,ire1samples]);
        fprintf('Now getting ire1 modules\n')
        [ire1_modules_up, ire1_modules_down]=getDiffModules(cid,pseudobulkexp,ire1samples,wtsamples);
        fprintf('Now getting ire1 modules\n')
        [ire1_modules_up2, ire1_modules_down2]=getDiffModules(cid,pseudobulkexp,ire1samples,[wtsamples,xbpsamples]);
        fprintf('Now getting ire1 modules\n')
        [mutant_modules_up, mutant_modules_down]=getDiffModules(cid,pseudobulkexp,[xbpsamples ire1samples],wtsamples);
        fprintf('Now getting ire1 modules\n')
	[xbp_modules_up3, xbp_modules_down3]=getDiffModules(cid,pseudobulkexp,xbpsamples,ire1samples);
         
	fid=fopen('merlin_module_annot_minsize10_beta_v2.txt','w');
        writeModules(ire1_modules_up,fid,'IRE1-KO-UP vs IRE1-WT_XBP1-WT');
        writeModules(ire1_modules_down,fid,'IRE1-KO-DOWN vs IRE1-WT_XBP1-WT');
        writeModules(xbp_modules_up,fid,'XBP1-KO-UP vs IRE1-WT_XBP1-WT');
        writeModules(xbp_modules_down,fid,'XBP1-KO-DOWN vs IRE1-WT_XBP1-WT');
        writeModules(mutant_modules_up,fid,'IRE1-KO_XBP1-KO-UP vs IRE1-WT_XBP1-WT');
        writeModules(mutant_modules_down,fid,'IRE1-KO_XBP1-KO-DOWN vs IRE1-WT_XBP1-WT');
	writeModules(xbp_modules_up2,fid,'XBP1-KO-UP vs IRE1-WT_XBP1-WT_IRE1-KO');
        writeModules(xbp_modules_down2,fid,'XBP1-KO-DOWN vs IRE1-WT_XBP1-WT_IRE1-KO');
        writeModules(ire1_modules_up2,fid,'IRE1-KO-UP vs IRE1-WT_XBP1-WT_XBP1-KO');
        writeModules(ire1_modules_down2,fid,'IRE1-KO-DOWN vs IRE1-WT_XBP1-WT_XBP1-KO');
        writeModules(xbp_modules_up3,fid,'XBP1-KO-UP vs IRE1-KO');
	writeModules(xbp_modules_down3,fid,'XBP1-KO-DOWN vs IRE1-KO');
	fclose(fid);
