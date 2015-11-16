############################ Third Version #################################

x <- as.numeric(Sys.getenv("SGE_TASK_ID"))

testEdit<-function(j){
	Data20121011<-read.table("Common_Heterozygous_with_info_filtered_4_regression_newVersion.txt_Editing_tested.txt_preprocessed", head=T, as.is=T, sep = "\t")
	Data20121011<-Data20121011[Data20121011$cRef!="-" & Data20121011$cVar!="-",]
	Data20121011<-Data20121011[Data20121011$cRef!="0" | Data20121011$cVar!="0",]
	VariantsID <- with(Data20121011, paste(Chromosome, Start, sep="."))
	Data20121011$VarID <- VariantsID
	Data20121011$cRef <- as.integer(Data20121011$cRef)
	Data20121011$cVar <- as.integer(Data20121011$cVar)
	Data20121011$Disease <- as.factor(Data20121011$Disease)
	Data20121011$Sample_name <- as.factor(Data20121011$Sample_name)

	NCov <-2 # what is this?????

	AllVariants <- unique(VariantsID)

	library(lme4)

	###Control how to split the file
	max <- 10000
	x<- seq_along(AllVariants)
	y<- x
	Split <- split(y, ceiling(x/max))

	Index <- Split[[j]]
	AllVariants.sub <- AllVariants[Index]
	###Control how to split the file ### end


	N <- length(AllVariants.sub)
	AllResults <- vector(length=N, mode="character")
	Mod0Warnings <- vector(length=N, mode="character")
	Mod1Warnings <- vector(length=N, mode="character")
	NA.test.type <- vector(length=N, mode="character")

	for (i in seq(N)){
		Indx<-Data20121011$VarID==AllVariants.sub[i] #chr15.100246942
		cat(i," ",AllVariants.sub[i],": ") 
		YDat<-data.frame(Ref=Data20121011$cRef[Indx], Var=Data20121011$cVar[Indx])
		
		if ((length(which(YDat$Ref>0))<2) || (length(which(YDat$Var>0))<2)) {# handle 0 count
			AllResults[i] <- "NA"
			NA.test.type[i] <- "0ReforVar"
			cat("0_NA","\n")
			next
		}

		XDat<-data.frame(Status=Data20121011$Disease[Indx], SampleID=as.factor(1:length(Indx[Indx])))
       
        if(length(XDat$Status) == length(XDat$Status[XDat$Status == "NOF"]) | length(XDat$Status) == length(XDat$Status[XDat$Status == "OA"])) {
                AllResults[i] <- "NA"
                NA.test.type[i] <- "1ConMissing"
                cat("MissCon_NA","\n")
                next
        }

	if(length( which(XDat$Status == "OA") ) < 6 & length( which(XDat$Status == "NOF") ) <= 3 ) {
            AllResults[i] <- "NA"
            NA.test.type[i] <- "TooFewSamples"
            cat("TooFewSamples_NA","\n")
            next
        }

		TotCounts<-rowSums(YDat)
		#if (Trace) cat(i," ",as.character(Vars[i]),": ")
		Y<-unlist(apply(YDat,1,function(x){rep(c(0,1),x)}))
		X<-as.data.frame(matrix(unlist(apply(cbind(TotCounts,XDat),1, function(x){rep(x[-1],as.integer(x[1]))})),byrow=T,ncol=NCov))
		WrkgDf<-cbind(Y,X)
		names(WrkgDf)<-c("Read",names(XDat)) 
		ConTab<-table(c(WrkgDf$Read,0,1),c(WrkgDf$Status,1,2) )-diag(c(1,1))
		TotReads<-sum(ConTab)
		DiagCounts<- ConTab[1,1]+ConTab[2,2]

		if (DiagCounts==0 || DiagCounts==TotReads) {
			P.val.fish <--fisher.test(ConTab)$p.value
			AllResults[i] <- P.val.fish
			NA.test.type[i] <- "fisher"
			cat("fisher", " ", P.val.fish,"\n")
        		next
		}

		H0.glm<-tryCatch(lmer(Read~(1|SampleID),data=WrkgDf),
                      error=function(e){return(NA)},
                      warning=function(w){
                                          Mod0Warnings[i] <- "Y"
                                          return(NA)}
                                          )
		HA.glm<-tryCatch(lmer(Read~Status+(1|SampleID),data=WrkgDf), 
                      error=function(e){return(NA)},
                      warning=function(w){
                                          Mod1Warnings[i] <- "Y"
                                          return(NA)}
                                          )
		if ((class(H0.glm)[1]=="mer") && (class(HA.glm)[1]=="mer")) {
			LogLike1<-logLik(HA.glm)[[1]]                      
			LogLike0<-logLik(H0.glm)[[1]]
			P.val=pchisq(2*(LogLike1-LogLike0),1,lower=F)
			NA.test.type[i] <- "glm"
			AllResults[i] <- P.val
			cat(P.val,"\n")
		} else {
			AllResults[i] <- "NA"
			NA.test.type[i] <- "Error_War"
			cat("Error_War"," ","NA","\n")
		}
			            
	} 
	table.name <- paste(paste("NewSlicingRestult/regression_array_test_slice", min(Index), max(Index), sep="_"),"txt", sep = ".")
	result <- data.frame(Variants_id=AllVariants.sub, P_value=AllResults, test_type = NA.test.type, Mod0_warnings = Mod0Warnings, Mod1_warnings = Mod1Warnings)
	write.table(result, file=table.name, row.name=F, col.name=F, quote=F,sep="\t")
}

testEdit(x)

#getRaw<-function(ToSelect,AllDat){
#	Indx<-!is.na(match(AllDat$VarID,ToSelect))
#	AllDat[Indx,]
#}


