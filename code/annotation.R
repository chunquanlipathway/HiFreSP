####################################################################
##get KO sub-pathway annotation
identifyGraph<-function(moleculeList,graphList,type="gene",background=getBackground(type),
   order="pvalue",decreasing=FALSE,locateOrg=TRUE,ignoreAmbiguousEnzyme=TRUE){
      if(typeof(moleculeList)!="character"){
	  print("warning: your moleculeList must be 'character' vector. Because the type of your current moleculeList is not correct, it has been conveted arbitrarily using the function as.character().")
	  moleculeList<-as.character(moleculeList)
	  }
      if(!exists("k2ri")) initializeK2ri()
	  graphList.length<-length(graphList)
	  if(graphList.length==0){
	     print("warning: there is no any graphs in graphList or these graphs are not available for pathway analyses.")
	  }	  
	  if(locateOrg==TRUE){
	     if(graphList.length>0){
	         gene2path<-get("gene2path",envir=k2ri)
	         org.path<-unique(as.character(gene2path[,2]))
		     org.path<-substring(org.path,nchar(org.path)-4,nchar(org.path))
		     graphList<-graphList[sapply(graphList,function(x) substring(x$number,0,5) %in% org.path)]
		 }
	  }
      annList<-list()
	  if(graphList.length>0){
      for(i in 1:length(graphList)){
            ann<-list(pathwayId=character(),pathwayName="not known",annMoleculeList=character(),annMoleculeNumber=0,
                      annBgMoleculeList=character(),annBgNumber=0,moleculeNumber=0,bgNumber=0,pvalue=1,fdr=1)

			ann$pathwayId<-paste("path:",graphList[[i]]$number,sep="")
            KOList<-unique(unlist(strsplit(unlist(strsplit(V(graphList[[i]])$names," ")),";")))
            #KOList<-unique(unlist(strsplit(V(graphList[[i]])$names,"[ ;]")))			
		if(type=="gene"||type=="gene_compound"){	
            if(graphList[[i]]$org=="ko"){		
              graphGeneList<-getGeneFromKO(KOList) 
			}
			else if(graphList[[i]]$org=="ec"){
			  graphGeneList<-getGeneFromEnzyme(KOList,ignoreAmbiguousEnzyme=ignoreAmbiguousEnzyme)
			}
			else{
			  org_idType<-unlist(strsplit(graphList[[i]]$org,";"))
			  if(org_idType[1]==getOrgAndIdType()[1]){
			     if(length(org_idType)==2){
				    if(org_idType[2]==getOrgAndIdType()[2]){
					    graphGeneList<-KOList
					}
					else{stop(paste("graph ",i,"  error: it is not ec, ko, or current org graph.",sep=""))}
				 }
				 else{
				     graphGeneList<-getGeneFromKGene(KOList)
				 }
			  }
			  else{stop(paste("graph ",i,"  error: it is not ec, ko, or current org graph.",sep=""))}
			}			
            
			
        }			
       if(type=="compound"||type=="gene_compound"){	
            graphCompoundList<-KOList[substring(KOList,0,5)=="cpd:C"] 
			graphCompoundList<-unique(substring(graphCompoundList,5))
	   }
	   if(type=="gene_compound"){
	        graphMoleculeList<-c(graphGeneList,graphCompoundList)
	   }else if(type=="gene"){
	        graphMoleculeList<-graphGeneList   
	   }else if(type=="compound"){
	        graphMoleculeList<-graphCompoundList  	   
	   }
       annotatedMoleculeList<-intersect(graphMoleculeList,moleculeList)
	   annotatedBackgroundList<-intersect(graphMoleculeList,background)	
            
            pathwayName<-graphList[[i]]$title
            if(length(pathwayName)!=0)
                ann$pathwayName<-pathwayName
            ann$annMoleculeList<-annotatedMoleculeList 
         
            ann$annMoleculeNumber<-length(annotatedMoleculeList)
			ann$annBgMoleculeList<-annotatedBackgroundList
            ann$annBgNumber<-length(annotatedBackgroundList)

            ann$moleculeNumber<-length(moleculeList)
            ann$bgNumber<-length(background)

            ann$pvalue<-1-phyper(ann$annMoleculeNumber-1,ann$annBgNumber,
                 ann$bgNumber-ann$annBgNumber,ann$moleculeNumber)
            
            annList[[i]]<-ann
      } 
	  }
	  if(length(annList)>0){
	     p_value<-sapply(annList,function(x) return(x$pvalue))
         #fdrtool.List<-fdrtool(p_value,statistic="pvalue",plot=FALSE,verbose=FALSE)
         	 
         #print(fdrtool.List$qval)
         #for(i in seq(annList)){
         #   annList[[i]]$qvalue<-fdrtool.List$qval[i]
		 #	annList[[i]]$lfdr<-fdrtool.List$lfdr[i]
         #}
		 print(p_value)
		 fdr.List<-p.adjust(p_value,"BH")
		 for(i in seq(annList)){
		     annList[[i]]$fdr<-fdr.List[i]
		 }
         #names(annList)<-sapply(graphList,function(x) x$number)
		 
         annList<-annList[sapply(annList,function(x) x$annMoleculeNumber>0)]
         #annList<-annList[order(sapply(annList,function(x) x[[order]]),decreasing=decreasing)]   
	  }
	  return(annList)	

}
#####################################################################
fdr.est<-function(p)
{
    m <- length(ind <- which(!is.na(p)))
    fdr <- rep(NA, length(p))
    stat <- cbind(1:length(p), p, fdr)
    stat[ind, 3] <- unlist(lapply(stat[ind, 2], function(x) {
        c <- length(which(stat[ind, 2] <= x))
        m * x/c
    }))
	##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    stat[ind, ] <- stat[ind, ][order(stat[, 2], decreasing = TRUE), 
        ]
    stat[ind, 3] <- cummin(stat[ind, 3])
    fdr <- stat[order(stat[, 1]), 3]
    fdr[which(fdr > 1)] <- 1
    return(fdr)
}
#####################################################################
printGraph<-function(ann,detail=FALSE){
	  if(detail==FALSE){
	  pathwayId<-sapply(ann,function(x) x$pathwayId)
      pathwayName<-sapply(ann,function(x) x$pathwayName)
      annMoleculeRatio<-sapply(ann,function(x) paste(x$annMoleculeNumber,x$moleculeNumber,sep="/"))
      annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
      pvalue<-sapply(ann,function(x) x$pvalue)
      #qvalue<-sapply(ann,function(x) x$qvalue)
	  #lfdr<-sapply(ann,function(x) x$lfdr)
	  fdr<-sapply(ann,function(x) x$fdr)
      #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annMoleculeRatio,
      #                       annBgRatio,pvalue,qvalue,lfdr))
      ann.data.frame<-data.frame(pathwayId=pathwayId,pathwayName=pathwayName,annMoleculeRatio=annMoleculeRatio,
	  annBgRatio=annBgRatio,pvalue=pvalue,fdr=fdr,stringsAsFactors=FALSE)							 
	  }
	  else{	 
      pathwayId<-sapply(ann,function(x) x$pathwayId)	  
	  pathwayName<-sapply(ann,function(x) x$pathwayName)
	  annMoleculeList<-sapply(ann, function(x){ paste(x$annMoleculeList,collapse=";") })
      annBgMoleculeList<-sapply(ann, function(x){ paste(x$annBgMoleculeList,collapse=";")})
	  annMoleculeRatio<-sapply(ann,function(x) paste(x$annMoleculeNumber,x$moleculeNumber,sep="/"))
      annBgRatio<-sapply(ann,function(x) paste(x$annBgNumber,x$bgNumber,sep="/"))
      pvalue<-sapply(ann,function(x) x$pvalue)
      #qvalue<-sapply(ann,function(x) x$qvalue)
	  #lfdr<-sapply(ann,function(x) x$lfdr)
	  fdr<-sapply(ann,function(x) x$fdr)
      #ann.data.frame<-as.data.frame(cbind(pathwayId,pathwayName,annMoleculeRatio,
      #                       annBgRatio,pvalue,qvalue,lfdr,annMoleculeList,annBgMoleculeList))
      ann.data.frame<-data.frame(pathwayId=pathwayId,pathwayName=pathwayName,annMoleculeRatio=annMoleculeRatio,
	  annBgRatio=annBgRatio,pvalue=pvalue,fdr=fdr,annMoleculeList=annMoleculeList,
	  annBgMoleculeList=annBgMoleculeList,stringsAsFactors=FALSE)								 
	  }
      return(ann.data.frame)
}