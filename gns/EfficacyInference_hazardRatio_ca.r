options(java.parameters="-DSIMULATION_THREAD_COUNT = 256")
options(java.mem = 64 * 1024)
library(REFSfs)
library(data.table)
library(ggplot2)




efficacyInference = function(
    training_data,
    ens,
    treatmentVariable, # "dosesReceived"
    subgroup_driver, # "lbNonDen.DENGUE_1_ACAMBIS_PUO_359_PRNT50_0_dosesReceived_preDose1"
    subgroup_outcome, # "totalHospitalDays_bothPhases"
    partitionIntervalLength=NULL # (age_granularity) interval size for riemann sum partition (I think...)
    numPartitionIntervals=20, # can set instead of partitionIntervalLength
){


    num_networks  =   fsNumNets(ens) 
    num_patients  =   nrow(training_data)#length(training_data[[1]])
    fixed_vars    =   fsFixedVars(ens)
    sampleSize    = 10

    ##*** CA: ## sequence (unique_age_vector) is used multiple times - set it once 
    subgroup_driver_range = range(training_data[[subgroup_driver]]) # age_range
    if(is.null(partitionIntervalLength) & !is.null(numPartitionIntervals)){
        partitionIntervalLength =  diff(subgroup_driver_range)/numPartitionIntervals # (age_granularity, interval size for riemann sum partition)
    }
    subgroup_driver_partition = seq(
          from=subgroup_driver_range[[1]], 
          to=(subgroup_driver_range[[2]]-partitionIntervalLength), 
          by=partitionIntervalLength
    )
    ##***

    #Let's generate the probability distribution of dm.AGE from the data
    cum_prob_greaterAge = c()
    ##*** CA: 
    # for (i in seq(from=age_range[1], to=(age_range[2]-age_granularity), by=age_granularity))
    for (i in subgroup_driver_partition)
    {
    # cum_prob_greaterAge = c(cum_prob_greaterAge, diff(ecdf(training_data[[subgroup_driver]])(c(i,i+age_granularity))))
    cum_prob_greaterAge = c(cum_prob_greaterAge, diff(ecdf(training_data[[subgroup_driver]])(c(i,i+partitionIntervalLength))))
    }
    ##***


    ##*** CA:
    # unique_age_vector          = seq(from=age_range[1], to=(age_range[2]-age_granularity), by=age_granularity)
    # cum_prob_greaterAge        = cbind(cum_prob_greaterAge, unique_age_vector)
    cum_prob_greaterAge        = cbind(cum_prob_greaterAge, subgroup_driver_partition)
    cum_prob_greaterAge        = as.data.frame(cum_prob_greaterAge)
    # names(cum_prob_greaterAge) = c('AgeProb',subgroup_driver)
    names(cum_prob_greaterAge) = c('Prob',subgroup_driver)
    ##***


    #The user specifies that Age ranges from 2 to 15, and they want a granularity of .5 age groups
    #Set up a dataframe with intervention values only.
    #We will merge this data with respective confounders later on, per network, to create baselineData for fsSimulateOverwrite
    intervention_baselineData  =  data.frame()
    # for (i in seq(from=age_range[1], to=(age_range[2]-age_granularity), by=age_granularity)) #4:30)
    for (i in subgroup_driver_partition) 
    {
        age_vector          = c(rep(i, 2*num_patients))
        treatment_vector    = c(rep(1, num_patients))
        no_treatment_vector = c(rep(0, num_patients))
        tr_notr             = c(treatment_vector, no_treatment_vector) 

        temp  = cbind(tr_notr, age_vector)
        intervention_baselineData = data.frame(rbind(intervention_baselineData, temp))
    }

    ## ** CA: changing name of dosesReceived variable
    # names(intervention_baselineData) = c("dosesReceived", subgroup_driver)
    names(intervention_baselineData) = c(treatmentVariable, subgroup_driver)
    ## **

    ##*** CA:
    # num_experiments                  = 2*length(seq(from=age_range[1], to=(age_range[2]-age_granularity), by=age_granularity))
    num_experiments                  = 2*length(subgroup_driver_partition)
    ## ***

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
        ## *** CA
        # dosesReceived_predictors = fsCausalVars(subens, c('dosesReceived'), cutoff=0, maxpath=1) 
        dosesReceived_predictors = fsCausalVars(subens, treatmentVariable, cutoff=0, maxpath=1) 
        ## *** 

        dm.Age.predictors        = fsCausalVars(subens, c(subgroup_driver), cutoff=0, maxpath=1) 
        #fsCausalVars(subens, c('subgroup_outcome')       , cutoff=0, maxpath=1) 
        condition_vars           = c(fixed_vars, dosesReceived_predictors, dm.Age.predictors)
        condition_vars           = unique(condition_vars)
        condition_baselineData   = training_data[, condition_vars]

        ## *** CA:
        # condition_baselineData$dosesReceived       = NULL
        condition_baselineData[[treatmentVariable]]       = NULL
        ## ***
        condition_baselineData[[subgroup_driver]]  = NULL

        condition_baselineData   = rbindlist(rep(list(condition_baselineData),num_experiments))
        baselineData             = cbind(condition_baselineData ,intervention_baselineData)
        baselineData             = as.data.frame(baselineData)
        simulation               = fsSimulateOverwrite(subens, c(subgroup_outcome), baselineData = baselineData , sampleSize=sampleSize, reportMeans = FALSE)
        simulation               = simulation$result[[1]]
        simulation$network       = NULL
        simulation$output        = NULL
        simulation$condition     = NULL
        #pmf1                     = simulation$pmf$"1"
        #pmf0                     = simulation$pmf$"0"
        #simulation$pmf           = NULL
        #simulation$pmf1          = pmf1
        #simulation$pmf0          = pmf0
        simulation               = cbind(simulation, sampleSize_ind)



        #Here we should compute P(vcd|T=1, Age>x) = SUM_i=x:max P(vcd|T=1, Age=i)*P(Age=i) and P(vcd|T=0, Age>x) for each x
        #Multiple those values by P(Age=x) from cum_prob_greaterAge
        simulation               = cbind(simulation, intervention_baselineData2)
        simulation$fixedDataRow  = NULL
        simulation               = as.data.table(simulation)

        ## *** CA:
        # simulation               = simulation[, lapply(.SD, mean), by = c('dosesReceived', subgroup_driver, 'sample')]
        simulation               = simulation[, lapply(.SD, mean), by = c(treatmentVariable, subgroup_driver, 'sample')]
        ## ***
        simulation               = as.data.frame(simulation)

        #get rid of alpha! use beta in place of lambda

        ## *** CA: for now, will rename beta as lambda, but should change this block (use setnames instead, etc)...
        if("beta" %in% colnames(simulation)){
            simulation$alpha = NULL
            names(simulation)        = c(treatmentVariable, subgroup_driver,'sample','lambda')
            drug_simulation          = simulation[simulation[[treatmentVariable]]==1,]
            no_drug_simulation       = simulation[simulation[[treatmentVariable]]==0,]
            names(drug_simulation)   = c(treatmentVariable,subgroup_driver,'sample','lambdaT1')
            names(no_drug_simulation)= c(treatmentVariable,subgroup_driver,'sample','lambdaT0')
            simulation               = merge(drug_simulation, no_drug_simulation, by=c(subgroup_driver,'sample'))
        }
        # simulation$dosesReceived.x = NULL
        # simulation$dosesReceived.y = NULL
        simulation[[paste0(treatmentVariable, ".x")]] = NULL
        simulation[[paste0(treatmentVariable, ".y")]] = NULL

        simulation                 = merge(simulation, cum_prob_greaterAge, by=subgroup_driver)

        #P(vcd=1|Age=i, T=y)*P(Age=i)
        ## *** CA:
        # simulation$Pvcd.X.PageT1   = simulation$lambdaT1 * simulation$AgeProb
        # simulation$Pvcd.X.PageT0   = simulation$lambdaT0 * simulation$AgeProb
        simulation$Pvcd.X.PageT1   = simulation$lambdaT1 * simulation$Prob
        simulation$Pvcd.X.PageT0   = simulation$lambdaT0 * simulation$Prob
        ## ***

        total = length(simulation[simulation$sample==1,][[1]])
        efficacy_network = data.frame()

        for (k in 1:sampleSize)
        {
            for(j in 1:total)
            {
                temp_pre         = simulation[simulation$sample==k,]
                temp             = temp_pre[j:total,]

                #For each x, SUM_i=x:max P(vcd=1|Age=i, T=y)*P(Age=i), grouped by treatment
                efficacy_calc    = (sum(temp$Pvcd.X.PageT1)/sum(temp$Pvcd.X.PageT0))
                ##*** CA
                # temp             = data.frame(efficacy_calc, unique_age_vector[j], i, k)
                temp             = data.frame(efficacy_calc, subgroup_driver_partition[j], i, k)
                ##***
                names(temp)      = c('efficacy', 'age_cutoff', 'network', 'sample')
                efficacy_network = data.frame(rbind(efficacy_network, temp))
            }
        }

        efficacy = data.frame(rbind(efficacy, efficacy_network))

        }

    # write.csv(efficacy, "~/Dengue/Dengue_Subgroup4.1.2/Dengue_SerotypeEfficacy_HospDays.csv",row.names = FALSE)
      return(efficacy)
}



runDir = "/projects/GNS-200415/models/sap_model1/cvRun_favor_txSwitch_v2/allObs"

efficacyResult = efficacyInference(
    training_data = fread(file.path(runDir, "train.csv"), data.table=F),
    ens = fsReadModel(file.path(runDir, "NetFragFile.txt.protobuf")),
    treatmentVariable = "treatmentArm.ACTARM", #replacing dosesReceived variable
    subgroup_driver   = "baselineClinicalLab.C__Reactive_Protein",
    subgroup_outcome  = "os.OSVAL",
    numPartitionIntervals = 5 ##
)





#stop here!

# efficacyc <- summarySE(efficacy, measurevar="efficacy", groupvars=c("age_cutoff"))

# #pd <- position_dodge(0.1) # move them .05 to the left and right

# file = '/home/dcunha/Dengue_SeroEfficacy_HospDays.pdf'
# pdf(file)
# ggplot(efficacyc, aes(x=age_cutoff, y=efficacy)) + 
#     geom_errorbar(aes(ymin=efficacy-ci, ymax=efficacy+ci), width=.1) +
#     geom_point()
# dev.off()





# summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
#                       conf.interval=.95, .drop=TRUE) {
#     library(plyr)

#     # New version of length which can handle NA's: if na.rm==T, don't count them
#     length2 <- function (x, na.rm=FALSE) {
#         if (na.rm) sum(!is.na(x))
#         else       length(x)
#     }

#     # This does the summary. For each group's data frame, return a vector with
#     # N, mean, and sd
#     datac <- ddply(data, groupvars, .drop=.drop,
#       .fun = function(xx, col) {
#         c(N    = length2(xx[[col]], na.rm=na.rm),
#           mean = mean   (xx[[col]], na.rm=na.rm),
#           sd   = sd     (xx[[col]], na.rm=na.rm)
#         )
#       },
#       measurevar
#     )

#     # Rename the "mean" column    
#     datac <- rename(datac, c("mean" = measurevar))

#     datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

#     # Confidence interval multiplier for standard error
#     # Calculate t-statistic for confidence interval: 
#     # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
#     ciMult <- qt(conf.interval/2 + .5, datac$N-1)
#     datac$ci <- datac$sd * ciMult

#     return(datac)
# }










