#Inference of Drivers of Differential Treatment Response
options(java.parameters="-DSIMULATION_THREAD_COUNT = 64")
options(java.mem = 50 * 1024)
library(REFSfs, lib.loc="/gnshealthcare/software/refs/fs/3.1.1/builds/30/")
library(plyr)
library(data.table)
#library(BMS)
library(ggplot2)
library(iPredictRE)




subpopulationInferenceEngine = function(directory, training_file, ensemble_file, drug_var, drug_levels, outcome_var){
#Protobuf location

## CA directory where run will live
## adding outcome var to the name because might run this for multiple outcomes
## name might be too long...
txEffect_directory = file.path(directory, paste0("iPredictRuns_txEffect_", outcome_var))


ens = fsReadModel(ensemble_file)
#BaselineData location
data = read.csv(training_file)
#Number of networks
number_networks = fsNumNets(ens)

#Experiment Data
experimentData = data.frame(experiment = as.character(c(drug_var)),
							           condition = as.character(c("no drug","drug")),
							           variable=as.character(c(drug_var)),
							           value=drug_levels)


#Run treatment intervention simulation for each network
outcome <- fsSimulateOverwrite(ens, c(outcome_var), baselineData=data, reportMeans=FALSE, sampleSize = 20, 
								experimentData=experimentData, conditionNames=as.character(c("no drug","drug"))) 
outcome = outcome$result[[1]]




#Compute the causal effect for each patient for each network
#############################################################################################################
#This part of the code needs accomodate different output variable distributions.
#Linear output should use the output column
#For discrete, we may want to use the pmf column when computing causal effect
#For survival, we may want to use the hazard(which will also be in the pmf column)
#
#Adjustment should be done when aggregating, and when computing the causal effect thereafter
#############################################################################################################
#avg_patient_network <- aggregate (output ~ condition + fixedDataRow + network, FUN = mean, data=outcome)
if("pmf" %in% colnames(outcome))
{

## CA rewriting to replace 0s with epsilon
epsilon = 1e-8
outcome = data.table(cbind(outcome[,1:4], outcome$pmf))
outcome[`0`==0, `0`:=epsilon]
outcome[,output:=`1`/`0`]# odds
# outcome[`0`==epsilon,.N,by=condition]
####
# pmf0 = outcome$pmf$"0"
# pmf1 = outcome$pmf$"1"
# odds = pmf1/pmf0
# outcome = cbind(outcome[,c("condition", "fixedDataRow", "network")], odds)
# names(outcome)[names(outcome) == 'odds'] <- 'output'
}



if("beta" %in% colnames(outcome))
{
outcome=data.table(outcome)
outcome[,output:=beta]
}

# outcome <- as.data.table(outcome)

avg_patient_network <- outcome[, lapply(.SD,mean), 
                           by = list(condition, fixedDataRow, network)]
avg_patient_network = as.data.frame(avg_patient_network)

drug_patient_network     = avg_patient_network[avg_patient_network$condition=='drug',]
no_drug_patient_network  = avg_patient_network[avg_patient_network$condition=='no drug',]
drug_patient_network$causalEffect  = drug_patient_network$output/no_drug_patient_network$output
drug_patient_network$causalEffect  = log(drug_patient_network$causalEffect)
ce_patient_network = drug_patient_network[,c("fixedDataRow", "network", "causalEffect")]


#Initialize per-network results data.frame
posterior_network_results = data.frame(causalPredictor=numeric(0),PosteriorInclusionProb=numeric(0),PosteriorExpectedEffect=numeric(0),Network=numeric(0))


data[] <- lapply(data, function(x) as.numeric(x))

#drops <- c("USUBJID","gtrends","weather.min_TemperatureF","weather.min_Humidity","weather.max_TemperatureF","weather.max_Humidity","weather.max_PrecipitationIn","worldBank.wbIndicator40_year_1","worldBank.wbIndicator72_year_1","worldBank.wbIndicator83_year_1","worldBank.wbIndicator26_year_1","worldBank.wbIndicator12_year_1","worldBank.wbIndicator55_year_1","worldBank.wbIndicator3_year_1","dm.COUNTRY.MYS","dm.COUNTRY.PHL","dm.COUNTRY.IND","dm.COUNTRY.VNM","dm.COUNTRY.THA")
#data = data[ , !(names(data) %in% drops)]




######
## CA adding directory for iPredict runs
	dir.create(file.path(txEffect_directory, "data"), recursive=T)
######

#Loop through each network, and regress vars that interact with the treatment pathway on the treatment effect.

## CA:
iPredictRuns = c()## save indices
##
for (i in 1:number_networks){
	temp_ens          = fsSubsetEnsemble(ens, i)
	causalVars        = fsCausalVars(temp_ens, c(outcome_var), cutoff = 0, maxpath = -1)
	downstreamVars    = fsDownstreamVars(temp_ens, c(drug_var), cutoff = 0, maxpath = -1)
	downstreamVars    = c(downstreamVars, drug_var)
	causal_predictors = causalVars[! causalVars %in% downstreamVars]
## CA commenting this out (drops not defined)
	# causal_predictors = causal_predictors[!causal_predictors%in%drops]
	## variable with & in name causes issues, remove for now 
	causal_predictors = causal_predictors[!causal_predictors%in%"baselineMedication.HERBALHOMEOPATHIC&_DIETARY_SUPPLEMENTS"]

##

## CA: moving if statment up - only creating treatmentEffect_data after if statment
	if (length(causal_predictors)==0| var(ce_patient_network[ce_patient_network$network==i,]$causalEffect)==0.0) {
		#print("next")
		next
		}

	else {	

## CA:
iPredictRuns = c(iPredictRuns, i)
## CA:

	treatmentEffect_data = data.frame(data[,causal_predictors])
	names(treatmentEffect_data) = c(causal_predictors)

	treatmentEffect_data$treatmentEffect = ce_patient_network[ce_patient_network$network==i,]$causalEffect
	column_order                         = c("treatmentEffect", causal_predictors)
	treatmentEffect_data                 = treatmentEffect_data[,column_order]
## CA - changing path to tx data
	write.csv(treatmentEffect_data, paste0(txEffect_directory,"/data/treatmentEffect_data",i,".csv"),row.names = FALSE)

##	

	
##
		iPredictRE.instance <- iPredictREInstance(
		 d.train=paste0(txEffect_directory, "/data/treatmentEffect_data",i,".csv"),
		 args=list('endpoints'="treatmentEffect",
		 'max-parameters'='3',
		 'max-terms'='2',
		 'max-degree'='2',
		 'maxent'='FALSE',
		 'separate-maxent'='0',
		 'ensemble-size'='1',
		 'annealing-schedule'='1,1,1', 
		 'moves'='1000',
		 'overlap'='.95',
		 'required-variables'='',
		 'generator'='beta',
		 'prior'='bic',
		 'verbosity'='3' 
		 ),

		 job.path=txEffect_directory,## CA - changing job path
		 job.name=paste0('iPredictRunDir',i), qsub=TRUE)

		iPredictRE.setup(iPredictRE.instance)
		run = iPredictRE.run(iPredictRE.instance,mem_requested=2,npc=1)
	}
}


## CA: adding system sleep for now (instead of checking that iPredict runs are done)
Sys.sleep(60)
##

### stop here!


## CA using r wrapper instead of bash. actually, iPredictRE.extractReport not writing csv for some reason...
lapply(iPredictRuns, function(runIndex){
	# load(file.path(directory, paste0("iPredictRuns_txEffect/iPredictRunDir",runIndex), "iPredictRE_Instance.RData"))
	# iPredictRE.extractReport(iPredictRE.instance, samples=30)## choosing samples arbitrarily
dir = file.path(txEffect_directory, paste0("iPredictRunDir",runIndex))
setwd(dir)
chkptFile = grep("chkp", list.files(dir), value=T)
chkptFile = chkptFile[nchar(chkptFile)==min(nchar(chkptFile))]# i think?
system(paste0("extractReportData -A ", chkptFile))
})

# for (i in 1:number_networks){
	#use extract report data to get mean effect size and inclusion probability
	# extractReportData -A  ipredictre.chkpt.p12783
	# command = './extractReportData_Automate.sh'
	# system(paste0(
	# 		'
	# 		#!/bin/bash; 
	# 		for i in {1..128}
	# 		do
	# 		cd /projects/GNS-200415/models/oakAndPop_model1/cvRun_npass100_Tstep5_128networks_10fold/allObs/iPredictRuns_txEffect/iPredictRunDir$i; 
	# 		pattern="ipredictre.chkpt*?"; file=( $pattern ); 
	# 		extractReportData -A  --checkpoint-file "${file[0]}"
	# 		done
	# 		'
	# 	))


# }


#Initialize per-network results data.frame
posterior_network_results = data.frame(causalPredictor=numeric(0),PosteriorInclusionProb=numeric(0),PosteriorExpectedEffect=numeric(0),Network=numeric(0))

## CA:
# for (i in 1:number_networks){
for (i in iPredictRuns){
###
	#directory      = '/home/dcunha/Dengue/Dengue_Subgroup3/'


	#For effect size, take the average of the "estimate" column grouped by the "term" column:
	#p15260_p_values.csv 
	#Place the average estimate per iPredict run into a file of 1024 average effect estimates
	#Plot histogram of effect size.
## CA: changing path
	# iPredictRuns_path = paste0(directory, 'iPredictRuns/iPredictRunDir',i)
	iPredictRuns_path = paste0(txEffect_directory, '/iPredictRunDir',i)
##

	command = paste0('cd ',iPredictRuns_path,';pattern=\"*p_value_summary.csv\"; echo $pattern;')
	file_name = system(command, intern=TRUE)
	effect_file_name = unlist(strsplit(file_name, split=' '))[1]
# grep("p_value_summary", list.files(iPredictRuns_path), value=T)
####


	effect_path = paste0(iPredictRuns_path, "/",effect_file_name)

	effect_file = read.csv(effect_path)

	#For inclusion probability, union all of the edge frequency files, then plot distribution.
	temp = data.frame(effect_file$term, effect_file$nnw, effect_file$mean.R, i)
	names(temp) = c("causalPredictor", "PosteriorInclusionProb", "PosteriorExpectedEffect","Network")
	posterior_network_results = data.frame(rbind(posterior_network_results, temp))

	#remove intercept
	posterior_network_results = posterior_network_results[!posterior_network_results$causalPredictor=='1',]
}

# # pdf(file = '/projects/GNS-200143/2stepcost_predictedHistogram.pdf')
# ggplot(posterior_network_results$PosteriorExpectedEffect, aes(vals, fill = measure)) + 
#     geom_histogram(alpha = 0.5, position = 'identity', binwidth = .1) +
#     labs(x = 'Effect of Age on Log(OddsRatio)', y = 'Frequency', 
#          title = 'Effect of Age on Log(OddsRatio) Histogram')
# dev.off()




#Produce E[PosteriorInclusionProb|D] for each variable 
#by summing the PosteriorInclusionProb for each variable, 
#and dividing by number of networks

#compute expected values; credible intervals will be computed below
posterior_results = aggregate(cbind(PosteriorInclusionProb, PosteriorExpectedEffect) ~ causalPredictor, FUN = sum, data=posterior_network_results)
posterior_results$PosteriorExpectedEffect = posterior_results$PosteriorExpectedEffect/as.numeric(number_networks)
posterior_results$PosteriorInclusionProb  = posterior_results$PosteriorInclusionProb/as.numeric(number_networks)
posterior_results                         = posterior_results[order(-posterior_results$PosteriorInclusionProb),]


#generate credible intervals
#Credible interval is mu +- sd/sqrt(num_nets)
#But there are only observations for when the variable is present.
#So the sd has to be adjusted by the number of zeros.
credible_intervals = data.frame(FifthPercentile_InclusionProb=numeric(0),NinetyfifthPercentile_InclusionProb=numeric(0),FifthPercentile_ExpectedEffect=numeric(0),NinetyfifthPercentile_ExpectedEffect=numeric(0))

for (i in posterior_results$causalPredictor){
	temp1 = posterior_network_results[posterior_network_results$causalPredictor==i,]
	temp1$causalPredictor=NULL
	temp1$Network=NULL

	number_zeros = number_networks - length(temp[,1])
	if(number_zeros>0.0){
		for (j in 1:number_zeros){
			temp1 = data.frame(rbind(temp1, c(0,0)))
		}
	}
	sd_IP = sd(temp1$PosteriorInclusionProb)
	sd_EE = sd(temp1$PosteriorExpectedEffect)

	se_IP = 2.0*sd_IP
	se_EE = 2.0*sd_EE

	temp2   = posterior_results[posterior_results$causalPredictor==i,]
	FP_IP  = temp2$PosteriorInclusionProb  - se_IP
	NFP_IP = temp2$PosteriorInclusionProb  + se_IP
	FP_EE  = temp2$PosteriorExpectedEffect - se_EE
	NFP_EE = temp2$PosteriorExpectedEffect + se_EE
	
	temp3 = data.frame(cbind(FP_IP,NFP_IP, FP_EE, NFP_EE))
	names(temp3) = c("FifthPercentile_InclusionProb", "NinetyfifthPercentile_InclusionProb", "FifthPercentile_ExpectedEffect","NinetyfifthPercentile_ExpectedEffect")
	credible_intervals = data.frame(rbind(credible_intervals, temp3))
	
}

#combine credible intervals with expected values
posterior_results = data.frame(cbind(posterior_results,credible_intervals))

## CA:
write.csv(posterior_results, file=file.path(txEffect_directory, "posterior_results.csv"), row.names=F)

return(posterior_results)
}




projectDir = "/projects/GNS-200415/models/" 
models = c("sap_model1", "sap_model2", "sap_model3", "oakAndPop_model1", "oakAndPop_model2")
runNames = "cvRun_favor_txSwitch_v2"
runDirs = file.path(projectDir, models, runNames, "allObs")
outcomes = c("os.OSVAL", "pfs.PFSVAL", "or.BORVAL")
########################
i = 5
directory = runDirs[[i]]
training_file = file.path(directory, "train.csv")
ensemble_file = file.path(directory, "NetFragFile.txt.protobuf")
drug_var = "treatmentArm.ACTARM"
drug_levels = c(0,1)
########################

## ## 30 october - deleting previous results (redoing for credible interval part)
## system(paste0("rm -r ", file.path(directory, "iPredictRuns_txEffect_*")))


posterior_results = lapply(outcomes, function(outcome){
	return(subpopulationInferenceEngine(
			directory=directory, 
			training_file=training_file, 
			ensemble_file=ensemble_file, 
			drug_var=drug_var, 
			drug_levels=drug_levels, 
			outcome_var=outcome
	))
})


########################################################################
## once all are done:
source("/home/calonso/generalFunctions.R")
names(runDirs) = models
runDirs = runDirs[!models %in% "sap_model3"]

txEffectDirs = paste0("iPredictRuns_txEffect_", outcomes)

allResultPaths = file.path(runDirs, rep(txEffectDirs, length(runDirs)), "posterior_results.csv")
allResults = lapply(allResultPaths, fread)

names(allResults) = gsub("\\/", "", gsub(
	"cvRun_favor_txSwitch_v2|allObs|iPredictRuns_|projects|GNS-200415|models|posterior_results.csv",     
	"", allResultPaths))

allResults = rbindlist(allResults, idcol=T)
allResults[,c("model", "outcome"):=tstrsplit(.id, "txEffect_")]
allResults[,.id:=NULL]

########################
# allResults[,quantile(PosteriorInclusionProb)]
# allResults[PosteriorInclusionProb>.8,]

allResults[,absExpectedEffect:=abs(PosteriorExpectedEffect)]
setorder(allResults, outcome, -absExpectedEffect)
# setorder(allResults, outcome, -PosteriorInclusionProb)
setcolorder(allResults, c("model", "outcome", names(allResults)[c(1:7, 10)]))

shortResults = allResults[,.(model, outcome, causalPredictor, PosteriorExpectedEffect, absExpectedEffect)]

shortResults[,.N,by=causalPredictor][order(-N)][N>1,]

quantiles = allResults[,list(value=quantile(PosteriorExpectedEffect)), by=.(model, outcome)]
quantiles[,quantile:=seq(.N), by=.(model, outcome)]
quantiles = transposeData(quantiles, lhs=c("model", "outcome"), rhs="quantile")
setorder(quantiles, outcome, model)

# allResults[,quantile(FifthPercentile_InclusionProb)]
# allResults[,quantile(NinetyfifthPercentile_InclusionProb)]
# allResults[FifthPercentile_InclusionProb>0,]
# allResults[!(FifthPercentile_ExpectedEffect<0 & NinetyfifthPercentile_ExpectedEffect>0),]

write.csv(allResults, file="/home/calonso/subpopResults.csv", row.names=F)

# scp -r @starmaster:/home/calonso/subpopResults.csv /Users/calonso/"OneDrive - GNS Healthcare"/Roche/.
############################################################################################


