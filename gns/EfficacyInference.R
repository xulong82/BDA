options(java.parameters="-DSIMULATION_THREAD_COUNT = 128")
options(java.mem = 50 * 1024)
library(REFSfs)
library(data.table)
library(ggplot2)

training_data =   read.csv("/home/dcunha/Dengue/Dengue_Subgroup4/train_09142017.csv")
ens           =   fsReadModel("/home/dcunha/Dengue/Dengue_Subgroup4/NetFragFile.txt.protobuf")

num_networks  =   fsNumNets(ens) 
num_patients  =   length(training_data[[1]])
fixed_vars    =   fsFixedVars(ens)
age_range     =   c(2,15)
age_granularity   = .5
subgroup_driver   = "dm.AGE"
subgroup_outcome  = "vcd.VCD_general"
sampleSize        = 10

#Let's generate the probability distribution of dm.AGE from the data
cum_prob_greaterAge = c()
for (i in seq(from=age_range[1], to=(age_range[2]-age_granularity), by=age_granularity))
{
cum_prob_greaterAge = c(cum_prob_greaterAge,diff(ecdf(training_data$dm.AGE)(c(i,i+age_granularity))))
}

unique_age_vector          = seq(from=age_range[1], to=(age_range[2]-age_granularity), by=age_granularity)
cum_prob_greaterAge        = cbind(cum_prob_greaterAge, unique_age_vector)
cum_prob_greaterAge        = as.data.frame(cum_prob_greaterAge)
names(cum_prob_greaterAge) = c('AgeProb','dm.AGE')



#The user specifies that Age ranges from 2 to 15, and they want a granularity of .5 age groups
#Set up a dataframe with intervention values only.
#We will merge this data with respective confounders later on, per network, to create baselineData for fsSimulateOverwrite
intervention_baselineData  =  data.frame()
for (i in seq(from=age_range[1], to=(age_range[2]-age_granularity), by=age_granularity)) #4:30)
{
age_vector          = c(rep(i, 2*num_patients))
treatment_vector    = c(rep(1, num_patients))
no_treatment_vector = c(rep(0, num_patients))
tr_notr             = c(treatment_vector, no_treatment_vector) 

temp  = cbind(tr_notr, age_vector)
intervention_baselineData = data.frame(rbind(intervention_baselineData, temp))
}
names(intervention_baselineData) = c("dosesReceived", "dm.AGE")
num_experiments                      = 2*length(seq(from=age_range[1], to=(age_range[2]-age_granularity), by=age_granularity))

#Later we have to repeat our intervention data to match our sampleSize
intervention_baselineData2 = intervention_baselineData[rep(seq_len(nrow(intervention_baselineData)), each=sampleSize),] #rbindlist(rep(list(intervention_baselineData),sampleSize))

#Set up sample indicator for aggregation purposes later on
sampleSize_ind           = data.frame(matrix(ncol = 1, nrow = sampleSize))
sampleSize_ind$sample    =  c(1:sampleSize)
sampleSize_ind           = rbindlist(rep(list(sampleSize_ind),num_experiments*num_patients))
sampleSize_ind[,1]       = NULL

efficacy =  data.frame()
for (i in 1:num_networks)
{
print(i)
subens                   = fsSubsetEnsemble(ens, i)
dosesReceived_predictors = fsCausalVars(subens, c('dosesReceived'), cutoff=0, maxpath=1) 
dm.Age.predictors        = fsCausalVars(subens, c('dm.AGE')       , cutoff=0, maxpath=1) 
#fsCausalVars(subens, c('subgroup_outcome')       , cutoff=0, maxpath=1) 
condition_vars           = c(fixed_vars, dosesReceived_predictors, dm.Age.predictors)
condition_vars           = unique(condition_vars)
condition_baselineData   = training_data[, condition_vars]

condition_baselineData$dosesReceived = NULL
condition_baselineData$dm.AGE        = NULL

condition_baselineData   = rbindlist(rep(list(condition_baselineData),num_experiments))
baselineData             = cbind(condition_baselineData ,intervention_baselineData)
baselineData             = as.data.frame(baselineData)
simulation               = fsSimulateOverwrite(subens, c(subgroup_outcome), baselineData = baselineData , sampleSize=sampleSize, reportMeans = FALSE)
simulation               = simulation$result[[1]]
simulation$network       = NULL
simulation$output        = NULL
simulation$condition     = NULL
pmf1                     = simulation$pmf$"1"
pmf0                     = simulation$pmf$"0"
simulation$pmf           = NULL
simulation$pmf1          = pmf1
simulation$pmf0          = pmf0
simulation               = cbind(simulation, sampleSize_ind)



#Here we should compute P(vcd|T=1, Age>x) = SUM_i=x:max P(vcd|T=1, Age=i)*P(Age=i) and P(vcd|T=0, Age>x) for each x
#Multiple those values by P(Age=x) from cum_prob_greaterAge
simulation               = cbind(simulation, intervention_baselineData2)
simulation$fixedDataRow  = NULL
simulation               = as.data.table(simulation)
simulation               = simulation[, mean(pmf1), by = list(dosesReceived, dm.AGE, sample)]
simulation               = as.data.frame(simulation)
names(simulation)        = c('dosesReceived','dm.AGE','sample','pmf1')


drug_simulation          = simulation[simulation$dosesReceived==1,]
no_drug_simulation       = simulation[simulation$dosesReceived==0,]
names(drug_simulation)   = c('dosesReceived','dm.AGE','sample','pmfT1')
names(no_drug_simulation)= c('dosesReceived','dm.AGE','sample','pmfT0')
simulation               = merge(drug_simulation, no_drug_simulation, by=c('dm.AGE','sample'))

simulation$dosesReceived.x = NULL
simulation$dosesReceived.y = NULL
simulation                 = merge(simulation, cum_prob_greaterAge, by='dm.AGE')

#P(vcd=1|Age=i, T=y)*P(Age=i)
simulation$Pvcd.X.PageT1   = simulation$pmfT1 * simulation$AgeProb
simulation$Pvcd.X.PageT0   = simulation$pmfT0 * simulation$AgeProb




total = length(simulation[simulation$sample==1,][[1]])
efficacy_network = data.frame()

for (k in 1:sampleSize)
{
for(j in 1:total)
{
temp_pre         = simulation[simulation$sample==k,]
temp             = temp_pre[j:total,]

#For each x, SUM_i=x:max P(vcd=1|Age=i, T=y)*P(Age=i), grouped by treatment
efficacy_calc    = (sum(temp$Pvcd.X.PageT0)-sum(temp$Pvcd.X.PageT1)) / sum(temp$Pvcd.X.PageT0)
temp             = data.frame(efficacy_calc, unique_age_vector[j], i, k)
names(temp)      = c('efficacy', 'age_cutoff', 'network', 'sample')
efficacy_network = data.frame(rbind(efficacy_network, temp))
}
}





efficacy = data.frame(rbind(efficacy, efficacy_network))

}

write.csv(efficacy, "./Dengue_Efficacy_Network.csv",row.names = FALSE)



read.csv("./Dengue_Efficacy_Network.csv")
efficacyc <- summarySE(efficacy, measurevar="efficacy", groupvars=c("age_cutoff"))

pd <- position_dodge(0.1) # move them .05 to the left and right

file = './Dengue_Efficacy_noPoly.pdf'
pdf(file)
ggplot(efficacyc, aes(x=age_cutoff, y=efficacy)) + 
    geom_errorbar(aes(ymin=efficacy-ci, ymax=efficacy+ci), width=.1) +
    geom_point()
dev.off()





summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$sd * ciMult

    return(datac)
}
