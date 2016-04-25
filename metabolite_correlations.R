# "m" is peaklist with metabolites in columns and samples in rows

pval<- matrix(data=NA, nrow=ncol(m), ncol=ncol(m))
Sp_cor <- matrix(data=NA, nrow=ncol(m), ncol=ncol(m))
rownames(pval) <- colnames(m)
colnames(pval) <- colnames(m)

rownames(Sp_cor) <- colnames(m)
colnames(Sp_cor) <- colnames(m)


for(i in 1:ncol(m)){
	for(j in 1:ncol(m))
	{
	cor <- cor.test(m[,i],m[,j],method=c("spearman"),exact=FALSE)
	pval[i,j]<-cor$p.value
	Sp_cor[i,j]<-cor$estimate
	}
}	

write.table(pval, file="met_cor_spearman_pval.txt", sep="\t", quote=F, col.names=NA)
write.table(Sp_cor, file="met_cor_spearman_rho_correlation.txt", sep="\t", quote=F, col.names=NA)

