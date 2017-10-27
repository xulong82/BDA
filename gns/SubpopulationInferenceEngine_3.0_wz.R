# Inference of Drivers of Differential Treatment Response
options(java.parameters="-DSIMULATION_THREAD_COUNT = 16")
library(REFSfs)
library(plyr)
library(data.table)
library(BMS)

subpopulationInferenceEngine = function(directory, drug_var, drug_levels, outcome_var, effect_type){
  
#Protobuf location
ens = fsReadModel(paste(directory, "NetFragFile.txt.protobuf", sep=""))
#BaselineData location
data = read.csv(paste(directory, "train.csv", sep=""))


#Experiment Data
experimentData = data.frame(experiment = as.character(c(drug_var)),
                            condition = as.character(c("no drug","drug")),
                            variable=as.character(c(drug_var)),
                            value=as.numeric(drug_levels)) #c(-1,1)))


#Run treatment intervention simulation for each network
outcome = fsSimulateOverwrite(ens, c(outcome_var), baselineData=data, reportMeans=FALSE, sampleSize = 20, experimentData=experimentData, conditionNames=as.character(c("no drug","drug"))) 
outcome = outcome$result[[1]]

# Compute the causal effect for each patient for each network
avg_patient_network <- aggregate (. ~ condition + fixedDataRow + network, FUN = mean, data=outcome) # average across 20 parameter draws
drug_patient_network     = avg_patient_network[avg_patient_network$condition=='drug',]
no_drug_patient_network  = avg_patient_network[avg_patient_network$condition=='no drug',]

# calculate causal effect based on scale parameter \beta from Weibull model

## causal effect for different outcomes are different
# survial: effect_type use beta
# binary or multinomial: use pmf at the specified level
# continuous: use output

if(effect_type == "Survival"){
  drug_patient_network$causalEffect  = log(drug_patient_network$beta) - log(no_drug_patient_network$beta) 
} else if(grepl("pmf", effect_type)){
  lvl = grep("[0-9]+", effect_type)
  drug_patient_network$causalEffect  = drug_patient_network$pmf[, lvl] - no_drug_patient_network$pmf[, lvl]
} else{
  drug_patient_network$causalEffect  = drug_patient_network$output - no_drug_patient_network$output
}

ce_patient_network = drug_patient_network[,c("fixedDataRow", "network", "causalEffect")]


#Initialize results data.frame
posterior_network_results = data.frame(causalPredictor=numeric(0),PosteriorInclusionProb=numeric(0),PosteriorExpectedEffect=numeric(0),Network=numeric(0))
number_networks = fsNumNets(ens)

#Retrieve causal vars for each network
for (i in 1:number_networks){
	temp_ens          = fsSubsetEnsemble(ens, i)
	causalVars        = fsCausalVars(temp_ens, c(outcome_var), cutoff = 0, maxpath = -1)
	downstreamVars    = fsDownstreamVars(temp_ens, c(drug_var), cutoff = 0, maxpath = -1)
	downstreamVars    = c(downstreamVars, drug_var)
	causal_predictors = causalVars[! causalVars %in% downstreamVars]


	if (length(causal_predictors)==0) {
		#print("next")
		next
		} else if (length(causal_predictors)==1){
		treatmentEffect_data = data.frame(data[,causal_predictors])
		names(treatmentEffect_data) = c(causal_predictors)

		treatmentEffect_data$treatmentEffect = ce_patient_network[ce_patient_network$network==i,]$causalEffect
		column_order                         = c("treatmentEffect", causal_predictors)
		treatmentEffect_data                 = treatmentEffect_data[,column_order]

		model_var = lm(treatmentEffect_data$treatmentEffect ~ 1 + as.numeric(unlist(treatmentEffect_data[2])))
		bic_intercept = BIC(lm(treatmentEffect_data$treatmentEffect ~ 1))/2.0
		bic_var = BIC(lm(treatmentEffect_data$treatmentEffect ~ 1 + as.numeric(unlist(treatmentEffect_data[2]))  ))/2.0
		
		PIP = exp(-bic_var)/(exp(-bic_var)+exp(-bic_intercept))
		Post.Mean = PIP*coef(model_var)[2][[1]]
		temp = data.frame(cbind(column_order[2],PIP, Post.Mean, i))
		names(temp) = c("causalPredictor", "PosteriorInclusionProb", "PosteriorExpectedEffect","Network")

		posterior_network_results = data.frame(rbind(posterior_network_results, temp))
	} else {
		treatmentEffect_data                 = data[,causal_predictors]	

		treatmentEffect_data$treatmentEffect = ce_patient_network[ce_patient_network$network==i,]$causalEffect
		column_order                         = c("treatmentEffect", causal_predictors)
		treatmentEffect_data                 = treatmentEffect_data[,column_order]

		bayesian_model = bms(treatmentEffect ~ ., data=treatmentEffect_data, mprior = "uniform", g = "UIP", user.int = F, mcmc="enumerate")
		co = data.frame(coef(bayesian_model))

		temp = data.frame(cbind(rownames(co),co$PIP,co$Post.Mean, i))
		names(temp) = c("causalPredictor", "PosteriorInclusionProb", "PosteriorExpectedEffect","Network")

		posterior_network_results = data.frame(rbind(posterior_network_results, temp))
	}
}


#Produce E[PosteriorInclusionProb|D] for each variable 
##by summing the PosteriorInclusionProb for each variable, 
##and dividing by number of networks

posterior_network_results$PosteriorInclusionProb  = as.numeric(as.character(posterior_network_results$PosteriorInclusionProb))
posterior_network_results$PosteriorExpectedEffect = as.numeric(as.character(posterior_network_results$PosteriorExpectedEffect))

posterior_results = aggregate(cbind(PosteriorInclusionProb, PosteriorExpectedEffect) ~ causalPredictor, FUN = sum, data=posterior_network_results, na.action = NULL)

posterior_results$PosteriorExpectedEffect = posterior_results$PosteriorExpectedEffect/as.numeric(number_networks)
posterior_results$PosteriorInclusionProb  = posterior_results$PosteriorInclusionProb/as.numeric(number_networks)
posterior_results                         = posterior_results[order(-posterior_results$PosteriorInclusionProb),]

return(posterior_results)
}




