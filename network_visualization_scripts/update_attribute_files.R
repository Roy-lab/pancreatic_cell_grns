
#in example of module 3339 this script updates the attribue file to change the TF (circle) in column2 to TF2 (diamond) if gene is not an enriched regulators in module, TAR (rectangle) if it is a module gene and will notchange the TF if the gene is an enriched regulators in module 

att=read.table("cytoscape_inputs/in_WT_xbp1_att.txt", header=T)
subnet=read.table("module3339_subnetwork.txt")
modgenes=read.table("genes_in_module_3339.txt")
enriched_genes=read.table("enriched_regulators_per_module.txt")
enriched_genes=enriched_genes[enriched_genes$V2 == "3339", ]

all_tfs=unique(as.character(subnet$V1))
all_tar=unique(as.character(subnet$V2))
modgenes=unique(as.character(modgenes$V1))
enriched=unique(as.character(enriched_genes$V1))
diamonds=setdiff(all_tfs,enriched)
circles=intersect(all_tfs,enriched)
rectangle=setdiff(modgenes,all_tfs)
att$Type=as.character(att$Type)
df=att
for (i in 1:nrow(df)) {if (df[i, 1] %in% diamonds) {df[i, 2] <- "TF2"} else if (df[i, 1] %in% circles) {df[i, 2] <- "TF"} else {df[i, 2] <- "TAR"}}

write.table(df,"cytoscape_inputs/in_WT_xbp1_2120_2118.txt",sep="\t",row.names=FALSE, quote=FALSE)