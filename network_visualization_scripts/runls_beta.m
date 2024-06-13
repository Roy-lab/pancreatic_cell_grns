function runls()

mods= {'2398','3061','3146','3339','3353','3476','3663','3793','3897','3916','3966','3483','3542','3781'}
ds={'KO_ire1_6680_6681','KO_xbp1_2117_2119','WT_ire1_6683','WT_xbp1_2120_2118'}

for i=1:14
	for j=1:4
		fprintf('%s\t%s\n',mods{i},ds{j})
		a = importdata(sprintf('../MERLIN_on_scVI_corrected_reseq_data/network_visualization/beta/module_%s/matlab_inputs/%s.txt',mods{i},ds{j}));
		names = a.textdata(:,1);
		%expression matrix
		train = a.data;
		train = zscore(train')';
		
		adj = load(sprintf('../MERLIN_on_scVI_corrected_reseq_data/network_visualization/beta/module_%s/matlab_inputs/adj.%s.txt',mods{i},ds{j}));
		adj = sparse(adj(:,1),adj(:,2),adj(:,3),size(train,1),size(train,1));
		adj = full(adj);
		%make sure there is no self loop
		for k=1:size(adj,1)
			adj(k,k)=0;
		end
		runOneNet(sprintf('../MERLIN_on_scVI_corrected_reseq_data/network_visualization/beta/module_%s/matlab_outputs/reg.%s.txt',mods{i},ds{j}),names,train,adj);
	end
end

function runOneNet(outname,names,train,adj)

fid = fopen(outname,'w');

for i=1:size(train,1)
	tfcount = sum(adj(i,:));
	if tfcount == 0
		%fprintf('%s\n',names{i});
		continue;
	end
	tfids = adj(i,:)~=0;
	%tf expression
	xx = train(tfids,:)';
	%tg expression
	yy = train(i,:)';
	v  = xx\yy;

	tfids = find(tfids);
	for j=1:length(tfids)
		fprintf(fid,'%s\t%s\t%f\n',names{tfids(j)},names{i},v(j));
	end
end
fclose(fid);
