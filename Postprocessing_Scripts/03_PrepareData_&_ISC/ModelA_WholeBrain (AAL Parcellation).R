#!/usr/bin/env AFNI_Batch_R
##!/usr/bin/env afni_run_R
library("brms")
library("hypr")
options(contrasts = c("contr.sum", "contr.poly"))# deviation coding: -1/0/1 
options(mc.cores = parallel::detectCores())
require('cmdstanr')
set_cmdstan_path(path = "//wsl$/Ubuntu-20.04/home/vasco/.cmdstan/cmdstan-2.36.0")

#Parse dems options
#initialize with defaults
lop <- NULL
lop$chains <- 4 #number of Markov chains
lop$iterations <- 2000 #number of iterations per Markov chain. Increase if Rhat greater than 1.1
lop$model  <- 'GroupStr*VideoType' #model equation
lop$cVars  <- 'GroupStr,VideoType' #categorical variables
lop$qVars  <- 'Intercept' #quantitative variables
lop$WCP    <- 10 #Within chain parallelization
lop$Y      <- 'Y' #response variable
lop$distY  <- 'student' #distribution for the response variable
lop$ROI   <- 'ROI' #ROI column naming
lop$Subj1   <- 'Subj1'
lop$Subj2   <- 'Subj2'
lop$r2z     <- TRUE # Fisher transformation of response variable
lop$outFN  <- paste0(getwd(),'/ModelA') #path prefix for output
lop$dataTable <- read.table(paste0(getwd(),'/data_aal_avg_modelA.txt'), header=T) #main data table
lop$overwrite <- TRUE                         
                                               
#Change options list to MBA variable list , check output files
process.MBA.opts <- function (lop, verb = 0) {
   if(is.null(lop$outFN)) stop(sprintf("Output filename not specified!"))
   if(isFALSE(lop$overwrite) && (
            file.exists(paste0(lop$outFN,".txt")) ||
            file.exists(paste0(lop$outFN,".RData")) ||
            file.exists(paste0(lop$outFN,".pdf"))) ) {
         stop(sprintf("File %s exists! Try a different name.", lop$outFN))
         return(NULL)
   }      
   if(!is.null(lop$cVars[1])) lop$CV <- strsplit(lop$cVars, '\\,')[[1]]   #split categorical variables
   if(!is.na(lop$qVars[1])) lop$QV <- strsplit(lop$qVars, '\\,')[[1]]     #split quantitative variables

   return(lop)
}
lop <- process.MBA.opts(lop, verb = lop$verb)

#################################################################################
########################## Begin MBA main #######################################
#################################################################################

# standardize the names for Y, ROI and subject
names(lop$dataTable)[names(lop$dataTable)==lop$ROI] <- 'ROI'
names(lop$dataTable)[names(lop$dataTable)==lop$Y] <- 'Y'
names(lop$dataTable)[names(lop$dataTable)==lop$Subj1] <- 'Subj1'
names(lop$dataTable)[names(lop$dataTable)==lop$Subj2] <- 'Subj2'

# make sure ROI1, ROI2 and Subj are treated as factors
if(!is.factor(lop$dataTable$Subj1)) lop$dataTable$Subj1 <- as.factor(lop$dataTable$Subj1)
if(!is.factor(lop$dataTable$Subj2)) lop$dataTable$Subj2 <- as.factor(lop$dataTable$Subj2)
if(!is.factor(lop$dataTable$ROI)) lop$dataTable$ROI <- as.factor(lop$dataTable$ROI)

# verify variable types
if(lop$model==1) terms <- 1 else terms <- strsplit(lop$model, split='\\*|\\+')[[1]]  #split the model by +, to ge the number of terms
if(length(terms) > 1) {
   for(ii in 1:length(terms)) {
      if(!is.null(lop$cVars[1])) if(terms[ii] %in% strsplit(lop$cVars, '\\,')[[1]] & !is.factor(lop$dataTable[[terms[ii]]])) # declared categorical variable with quantitative levels
         lop$dataTable[[terms[ii]]] <- as.factor(lop$dataTable[[terms[ii]]])
      if(terms[ii] %in% strsplit(lop$qVars, '\\,')[[1]] & is.factor(lop$dataTable[[terms[ii]]])) # declared quantitative variable contains characters
         stop(sprintf('Column %s in the data table is declared as numerical, but contains characters!', terms[ii]))
   }
}

# Fisher transformation
fisher <- function(r) ifelse(abs(r) < .995, 0.5*(log(1+r)-log(1-r)), stop('Are you sure that you have correlation values so close to 1 or -1?'))
if(lop$r2z) lop$dataTable$Y <- fisher(lop$dataTable$Y)

###Constrasts###

#Make sum coding Group
HcSumGroup <- hypr(
  Group1 = (HC_HC + SCZA_SCZA + SCZB_SCZB) / 3 ~ 0,
  Group2 = SCZA_SCZA ~ HC_HC,
  Group3 = SCZB_SCZB ~ HC_HC,
  levels = c("HC_HC", "SCZA_SCZA", "SCZB_SCZB")
)
contrasts(lop$dataTable[["GroupStr"]])<-contr.hypothesis(HcSumGroup)

#Make sum coding Socialness
HcSumType <- hypr(
  Type1 = (neutral + socialpos + socialneg + erotic) / 4 ~ 0,
  Type2 = neutral ~ (socialpos + socialneg + erotic) / 3,
  Type3=erotic~socialneg,
  Type4=erotic~socialpos,
  levels = c('erotic','neutral','socialneg','socialpos')
)
contrasts(lop$dataTable[["VideoType"]])<-contr.hypothesis(HcSumType)

##### model formulation #####
set.seed(1234)
lop$dataTable$w <- 1
modelForm <- as.formula(paste('Y~', lop$model, '+(', lop$model, '| gr(ROI, dist = "student"))+(VideoType| Subj1:Subj2)+(VideoType|mm(Subj1, Subj2, weights = cbind(w,w), scale=FALSE))'))


#(1|mm(ROI:Subj1, ROI:Subj2, weights = cbind(w, w), scale=FALSE))')





##SetPriors##
priors <-c(prior(student_t(3, 0, 0.2), class = "Intercept"), #Intercept
           prior(student_t(3, 0, 0.2), class = "sd"), #StdDev
           prior(student_t(3, 0, 0.3),class="sigma"), #Residuals
           prior(student_t(3, 0, 0.1), class = "b")) #Effects

##################### MCMC ####################
#
set.seed(1234)
ptm <- proc.time()
fm <- brm(modelForm, data=lop$dataTable, family=lop$distY, 
         prior=priors,
         chains = lop$chains, iter=lop$iter, warmup=1000, 
         backend = "cmdstanr", threads = threading(lop$WCP), refresh=1, control = list(max_treedepth = 15), adapt_delta = 0.99)
print(format(Sys.time(), "%D %H:%M:%OS3"))
proc.time() - ptm

#
set.seed(1234)
ptm <- proc.time()
fm <- brm(modelForm, data=lop$dataTable, family=lop$distY, 
          prior=priors,
          chains = lop$chains, iter=lop$iter, warmup=1000, 
          backend = "cmdstanr", threads = threading(lop$WCP), refresh=100, sample_prior = 'only')
print(format(Sys.time(), "%D %H:%M:%OS3"))
proc.time() - ptm


set.seed(1234)
ptm <- proc.time()
vm <- brm(modelForm, data=lop$dataTable, family=lop$distY, prior=priors,
          iter=10000, backend = "cmdstanr", threads = threading(16), 
          algorithm="fullrank", output_samples = 4000)
print(format(Sys.time(), "%D %H:%M:%OS3"))
proc.time() - ptm


ptm <- proc.time()
vm <- brm(modelForm, data=lop$dataTablefiltered, family=lop$distY, prior=priors,
          backend = "cmdstanr", threads = threading(16), 
          algorithm="pathfinder", draws = 4000, chains=10,psis_resample=TRUE,
          max_lbfgs_iters=10000, history_size=10)
print(format(Sys.time(), "%D %H:%M:%OS3"))
proc.time() - ptm
print(format(Sys.time(), "%D %H:%M:%OS3"))
proc.time() - ptm


save.image(file=paste0(getwd(), "/FullStdNoSubjStdPrior.RData"))


random_effects <- ranef(fm4, summary = FALSE) 
roiLabels <- dimnames(random_effects$ROI)[[2]] #ROI labels
conditional_effects(fm4, effects='GroupStr', conditions=roiLabels)



fm <- brm(bf(modelForm, decomp = 'QR'), data=lop$dataTable, family='student',
          prior=priors,
          chains = lop$chains, iter=2000, warmup=1000,
          backend = 'cmdstanr', threads = threading(lop$WCP, force=TRUE), refresh=1,
          stan_model_args = list(stanc_options = list('O1')),
          normalize=FALSE)
