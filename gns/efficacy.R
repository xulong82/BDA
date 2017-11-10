options(java.parameters="-DSIMULATION_THREAD_COUNT = 128")
options(java.mem = 50 * 1024)
library(REFSfs)
library(data.table)
library(ggplot2)

load("~/roche34/otf2/alldata.rdt") # training data
training_data = refsdf.c1

ens = fsReadModel("/users/xwang/roche34/otf2/m1.new/NetFragFile.txt.protobuf") # network ensemble
edges = fsEdgeFrequencies(ens, freqThreshold = 0) # edges

num_networks  =   fsNumNets(ens) 
num_patients  =   length(training_data[[1]])
fixed_vars    =   fsFixedVars(ens)

var = "asl.clin.num_EXINV2"
table(training_data[var])

var_range     =   c(0, 6)
var_granularity   = 1
subgroup_driver   = var
subgroup_outcome  = "pfs_AVAL"
sampleSize        = 10

# probability distribution of the variable of interest from training data
cum_prob = c() # ? cpdf or pdf
for (i in seq(from=var_range[1], to=(var_range[2]-var_granularity), by=var_granularity))
{
cum_prob = c(cum_prob, diff(ecdf(training_data[, var])(c(i,i+var_granularity))))
}

unique_var_vector          = seq(from=var_range[1], to=(var_range[2]-var_granularity), by=var_granularity)
cum_prob_greater        = cbind(cum_prob, unique_var_vector)
cum_prob_greater        = as.data.frame(cum_prob_greater)
names(cum_prob_greater) = c('varProb', var)

# set up a dataframe with intervention values only
# we will merge this data with respective confounders later on, per network, to create baselineData for fsSimulateOverwrite
# ? per patient, per experiment condition, per biomarker level 
intervention_baselineData  =  data.frame()
for (i in seq(from=var_range[1], to=(var_range[2]-var_granularity), by=var_granularity)) #4:30)
{
var_vector          = c(rep(i, 2*num_patients))
treatment_vector    = c(rep(1, num_patients))
no_treatment_vector = c(rep(0, num_patients))
tr_notr             = c(treatment_vector, no_treatment_vector) 

temp  = cbind(tr_notr, var_vector)
intervention_baselineData = data.frame(rbind(intervention_baselineData, temp))
}
names(intervention_baselineData) = c("dosesReceived", var)
num_experiments                      = 2*length(seq(from=var_range[1], to=(var_range[2]-var_granularity), by=var_granularity))

# later we have to repeat our intervention data to match our sampleSize
intervention_baselineData2 = intervention_baselineData[rep(seq_len(nrow(intervention_baselineData)), each=sampleSize),] 
# rbindlist(rep(list(intervention_baselineData),sampleSize))

#Set up sample indicator for aggregation purposes later on
sampleSize_ind           = data.frame(matrix(ncol = 1, nrow = sampleSize))
sampleSize_ind$sample    =  c(1:sampleSize)
sampleSize_ind           = rbindlist(rep(list(sampleSize_ind),num_experiments*num_patients))
sampleSize_ind[,1]       = NULL

var_treat = "asl.treat.char_ARMCD"

efficacy =  data.frame()
for (i in 1:num_networks)
{
print(i)
subens                   = fsSubsetEnsemble(ens, i)
dosesReceived_predictors = fsCausalVars(subens, c(var_treat), cutoff=0, maxpath=1) 
dm.var.predictors        = fsCausalVars(subens, c(var)       , cutoff=0, maxpath=1) 
#fsCausalVars(subens, c('subgroup_outcome')       , cutoff=0, maxpath=1) 
condition_vars           = c(fixed_vars, dosesReceived_predictors, dm.var.predictors)
condition_vars           = unique(condition_vars)
condition_baselineData   = training_data[, condition_vars]

condition_baselineData$dosesReceived = NULL
condition_baselineData$dm.var        = NULL

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

# compute P(vcd|T=1, Age>x) = SUM_i=x:max P(vcd|T=1, Age=i)*P(Age=i) and P(vcd|T=0, Age>x) for each x
# multiple those values by P(Age=x) from cum_prob_greaterAge
simulation               = cbind(simulation, intervention_baselineData2)
simulation$fixedDataRow  = NULL
simulation               = as.data.table(simulation)
# simulation               = simulation[, mean(beta), by = list(dosesReceived, "asl.clin.num_EXINV2", sample)] # ???
# simulation               = simulation[, list(alpha = mean(alpha), beta = mean(beta)), by = "dosesReceived,asl.clin.num_EXINV2,sample"]
simulation               = as.data.frame(simulation)
names(simulation)        = c('dosesReceived',var,'sample','beta')

drug_simulation          = simulation[simulation$dosesReceived==1,]
no_drug_simulation       = simulation[simulation$dosesReceived==0,]
names(drug_simulation)   = c('dosesReceived',var,'sample','beta1')
names(no_drug_simulation)= c('dosesReceived',var,'sample','beta0')
simulation               = merge(drug_simulation, no_drug_simulation, by=c(var,'sample'))

simulation$dosesReceived.x = NULL
simulation$dosesReceived.y = NULL
simulation                 = merge(simulation, cum_prob, by=var)

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
