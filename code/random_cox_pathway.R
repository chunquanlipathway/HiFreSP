rm(list=ls())
##################reading the gene expression profile ########################################
genedata <- read.table("input/genedata_eg.csv",header=T,sep=",") 

#################setting the random_times ######################
random_times <- 100
###################################################################
sample_num <- length(genedata[,1])

library(survival)
library(iSubpathwayMiner);
idconvert <- read.table("input/convertable.csv",header=TRUE,sep=",")
idconvert<-na.omit(idconvert)

#######################################################################
#######convert all metabolic pathways to graphs with KOs and compounds as nodes.
g1<-getMetabolicKOKOUGraph();   
g2<-getNonMetabolicKOKOUGraph(); 
graphList<-c(g1,g2);

endvalue <- length(genedata[1,])  
trnsave <- matrix(c(1:sample_num),nc=1)
x <- c(1:sample_num)
mRNArbst<-c()
pid <- c()

################ bootstrap processing ############################
for(j in 1:random_times)
{

	y <- sample(x,sample_num,replace=T)
	bootstap_dataset_1 <- genedata[y,]
	trnsave <- cbind(trnsave,bootstap_dataset_1[,1])

#######################  cox univariable analysis  #########################

	result_mRNA <- data.frame(colnames(genedata),NA,NA) 
	colnames(result_mRNA) <- c('mRNA','coef','p')
	for(i in 4:endvalue)     
	{
   		res <- coxph(Surv(bootstap_dataset_1$st,bootstap_dataset_1$status)~as.numeric(bootstap_dataset_1[,i]))
   		result_mRNA[i,2] <- res$coef  
   		result_mRNA[i,3]<- summary(res)$coefficients[5]   
	}
	result_mRNA <- result_mRNA[4:length(result_mRNA[,1]),]  
	write.csv(result_mRNA,paste("output/coxresult/coxresult_",j,".csv",sep=''),row.names=F)
	
#######################  pathway enrichment  #########################

	diff_genes<-result_mRNA[result_mRNA$p<0.05,1]
	length(diff_genes)  #2708

	locations<-match(diff_genes,idconvert[,1]);  
	diff_geneid<-idconvert$EntrezGene.ID[locations];

	geneList <-as.character(diff_geneid)  

	#and identify significant subpathways
	ann<-identifyGraph(geneList,graphList,type="gene",background=getPrioBackground(type="gene"))
	result<-printGraph(ann,detail=TRUE)
	result_risk<-result[which(result[,"pvalue"]<0.05),]

	##write.table(result_risk,file=paste("output/enrichment_pathway/all/all_subpathways_random_sample",j,".txt",sep=''),row.names=FALSE,sep="\t")
	write.table(result_risk,file=paste("output/enrichment_pathway/signifi_pathways_random_sample",j,".txt",sep=''),row.names=FALSE,sep="\t")

	
	ff<- print(ann)
	endvalue2 <- length(result[,1])
	
	for(k in 1:endvalue2)
	{
		if(ff[[k]]$pvalue < 0.05)
		{
			mRNArbst<-c(mRNArbst,ff[[k]]$annMoleculeList)
			pid <- c(pid,ff[[k]]$pathwayId)
			#number <- c(number,ff[[k]]$annMoleculeNumber)
		}
	}

}

### save the bootstrap sample ID ###################
#write.csv(trnsave,"output/random_sam_save.csv",row.names=F)


#write.table(mRNArbst,file="output/mRNArbst.txt",row.names=FALSE,sep="\t")
#write.table(pid,file="output/pid.txt",row.names=FALSE,sep="\t")
#write.table(table(mRNArbst),file="output/mRNArbst_statis.txt",row.names=FALSE,sep="\t")
#write.table(table(pid),file="output/pid_statis.txt",row.names=FALSE,sep="\t")


################### significant pathway information£¨p<0.05£©#############
risk_all<- read.table(paste("output/enrichment_pathway/signifi_pathways_random_sample",1,".txt",sep=''),header=TRUE,sep="\t")


for(i in 2:random_times)
{
	risk_p <- read.table(paste("output/enrichment_pathway/signifi_pathways_random_sample",i,".txt",sep=''),header=TRUE,sep="\t")	
	risk_all <- rbind(risk_all,risk_p) 

}
dim(risk_all)

#write.csv(risk_all,"output/risk_all_risk.csv",row.names=T)
pid_statis <- table(as.character(risk_all$pathwayId))
write.table(pid_statis,file="output/pid_statis.txt",row.names=FALSE,sep="\t")


pid_names<- read.table("output/pid_statis.txt",header=TRUE,sep="\t")
#pid_names<-pid_statis
endvalue <- length(pid_names[,1])
for(i in 1:endvalue)
{
	x1 <- as.character(risk_all[which(as.character(risk_all$pathwayId)==as.character(pid_names[i,1])),7])
	z <- c()
	for(j in 1:pid_names[i,2])
	{
		y <- unlist(strsplit(x1[j],";"))
		z <- c(z,y)		
	}
	pna <- substr(as.character(pid_names[i,1]),6,10)
	write.table(table(z),paste("output/pmRNA/statis_path",pna,".txt",sep=''),row.names=FALSE,sep="\t")
}

