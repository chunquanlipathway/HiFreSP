rm(list=ls())

###### dependent functions #####
source("annotation.R")

##################reading the gene expression profile ############

train_set <- read.table("input/genedata_eg.csv",header=TRUE,sep=",") 

#################setting the random_times ######################
random_times <- 100

#### set the threshold of HFP score  ####################

pid_names<- read.table("output/pid_statis.txt",header=TRUE,sep="\t")
pid_names_decr <- pid_names[order(pid_names[,2],decreasing=TRUE),]

T <- 0.5*max(pid_names$Freq)


##############################################################################
library(iSubpathwayMiner);
library(survival)

idconvert <- read.table("input/convertable.csv",header=TRUE,sep=",")
idconvert<-na.omit(idconvert)



g1<-getMetabolicKOKOUGraph();   
g2<-getNonMetabolicKOKOUGraph(); 
graphList<-c(g1,g2);

pid_names_2 <- pid_names[which(pid_names$Freq>T),] 



y <- matrix(c("pathwayid","subpathwayid","r_mRNA","ptrain_KM"),nc=4)

for(k in 1:length(pid_names_2[,1]))
{
	pna <- substr(as.character(pid_names_decr[k,1]),6,10)
	robs_mRNA <- read.table(paste("output/pmRNA/statis_path",pna,".txt",sep=''),header=TRUE,sep="\t")
	robs_mRNA <- robs_mRNA[robs_mRNA$Freq > 0.05*random_times,1]
	
#################  mining sub-pathway #############################

	geneList <- as.character(robs_mRNA)
	#####################################################################
    #############locate subpathways  with n and  s ########################
	#####################################################################
	subGraphList<-getLocSubGraph(geneList,graphList[pna],type="gene",n=1,s=5)
	
	if(length(subGraphList)>0)
	{
		ann<-identifyGraph(geneList,subGraphList,type="gene",background=getPrioBackground(type="gene"))
		result<-printGraph(ann,detail=TRUE)
		result_risk <- result
		write.table(result,file=paste("output/subpathway/subpathways_path",pna,".txt",sep=''),row.names=FALSE,sep="\t")

		for(q in 1:length(result_risk[,1]))
		{
			if(pna==substr(as.character(result_risk[q,1]),6,10))
			{
				r_mRNA <- unlist(strsplit(result_risk[q,7],";"))
				gene_candi <- idconvert[idconvert$EntrezGene.ID %in% r_mRNA,1]
				
                ################ risk score in the dataset ###################
				if(length(gene_candi)>0)
				{
					esca <- train_set
					#cols <-substring(colnames(esca)[11:12947],2, ); 
					#colnames(esca)[11:12947] <- cols
					esca_t <- esca[,colnames(esca) %in% as.character(gene_candi)]
					esca_t <- cbind(esca[,1:3],esca_t)
					if(length(gene_candi)==1)
					{
						colnames(esca_t)[4] <- gene_candi
					}
                    ########## cox univariate analysis ############
					
					result <- data.frame(colnames(esca_t),NA,NA) # prepare the result file 
					colnames(result) <- c('mRNA','coef','p')

					endvalue <- length(esca_t[1,])  #26543
					for(i in 4:endvalue)     
					{
						res <- coxph(Surv(esca_t$st,esca_t$status)~as.numeric(esca_t[,i]))
						result[i,2] <- res$coef  # beta 
						result[i,3]<- summary(res)$coefficients[5]   # p value
					}
					result <- result[4:length(result[,1]),]  

					rs_exp <- esca_t
					rs_all <- c()
					for(i in 1:length(esca_t[,1]))
					{
						for(j in 4:endvalue)
						{
							coef <- result[which(result[,1]==colnames(rs_exp)[j]),2]
							rs_exp[i,j] <- coef*rs_exp[i,j]
						}
						rs <-sum(rs_exp[i,4:endvalue])
						rs_all <- c(rs_all,rs)
					}


                    ############### KM result ############################
					re<-c()

					Y2 <- Surv(esca_t$st,esca_t$status)
					gene1<-rs_all
					
					gene2<-replace(gene1,gene1>=mean(gene1),1)
					gene2<-replace(gene2,gene1<mean(gene1),0)

					if((length(gene2[gene2==0])!=0)&(length(gene2[gene2==1])!=0))
					{
						survfit(Y2~gene2)
						sdf <- survdiff(Y2~gene2,rho=0)     #With rho = 0 this is the log-rank or Mantel-Haenszel test, and with rho = 1 it is equivalent to the Peto & Peto modification of the Gehan-Wilcoxon test. 
						p2 <- pchisq(sdf$chisq, df=1, lower=FALSE)
						re<-c(re,p2)
					}else{p2 <- 2}
					ptrain <-p2

					pdf(file=paste("output/pdf/train_path",pna,substr(as.character(result_risk[q,1]),11,13),".pdf",sep=''))#### plot KM curve
					plot(survfit(Y2~gene2), lty = c("solid", "solid"), col = c("blue", "red"),lwd=2,main=as.character(result_risk[q,1]), xlab="Time(months)",ylab="Overall survival")
					text(20, 0.05, paste("p=",as.character(round(p2,20))))
					legend("topright", c("low", "high"), col = c("blue", "red"), lty = c("solid", "solid")) 
					dev.off()

					###cutoff <- mean(gene1)

					y1 <- cbind(as.character(pid_names_decr[k,1]),as.character(result_risk[q,1]),result_risk[q,7],ptrain)

					y <- rbind(y,y1)
				}


			}
		}
	}		

}

write.table(y,"output/result_subpathway_n1s5.txt",row.names=FALSE,sep="\t")

