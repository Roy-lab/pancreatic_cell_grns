
clusters={'beta'};
sampleloc{1}=[446,300,796,548,718,1986,2775];
samplenames={'KO_ire1_6680'; 'KO_ire1_6681'; 'WT_ire1_6683';'KO_xbp1_2117';'KO_xbp1_2119';'WT_xbp1_2120';'WT_xbp1_2118'};
isletmarkernames={'Ins1';'Ins2';'Sst';'Gcg';'Ucn3';'Pdx1';'Slc2a2'};
for ind=1
cluster=clusters{ind};
DIRNAME=sprintf('expression_matrices');
CLUSTERDIR=sprintf('%s/modules',DIRNAME);
CLUSTERDIR_REG=sprintf('%s/list_of_terms',DIRNAME);
EXPDIR='beta_cells/'
EXPFILE=sprintf('%s/beta_cells_expression_matrix.txt',EXPDIR);
EXPFILE
alldata=importdata(EXPFILE);
gnames=alldata.textdata(2:end,1);
alldata=alldata.data;
figure;
t=0.4   
        CASSIGN=sprintf('%s/module_%s.%.1f_geneset.txt',CLUSTERDIR,cluster,t);
        cid=importdata(CASSIGN);
        markers=cid.textdata;

        CASSIGN_reg=sprintf('%s/list_%s.txt',CLUSTERDIR_REG,cluster);
        cid_reg=importdata(CASSIGN_reg);
        markers_reg=cid_reg.textdata; % gene names of the cluster assignment file
        [gid,actualid]=getGeneIDs(gnames,markers);
        expdata=alldata(gid,:);
        [gid_reg,actualid_reg]=getGeneIDs(gnames,markers_reg);
        expdata_reg=alldata(gid_reg,:);
        fprintf('Found %d genes\n',length(gid));
        outfname=sprintf('results/beta_cells%s_t%.1f',cluster,t);
        mkdir(outfname);
        [gid_isletmarker,actualid_islet]=getGeneIDs(gnames,isletmarkernames);
        expdata_islets=alldata(gid_isletmarker,:);
        figure;
        csize=showClusterWithReg_All_IsletMarker(cid,cid_reg,expdata,expdata_reg,4,sampleloc{ind},samplenames,outfname,expdata_islets,isletmarkernames);
end 
