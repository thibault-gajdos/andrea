rm(list=ls(all=TRUE))  ## efface les données
source('~/thib/projects/tools/R_lib.r')
setwd('~/thib/projects/andrea_stan/')


## * Prepare data
subject = 1
data <- read_csv('ares_expt1_bhv.csv')
data <- data %>% filter(subject == subject) %>%
    mutate(SOAactual = SOAactual * 10)

T <- nrow(data)
data_list <- list(
    T        = T,
    std_dev  = sd(data$SOAactual),
    delay    = data$SOAactual,
    choice   = data$synchrony
    )

## * fit
modelFile <- '~/thib/projects/andrea_stan/asym_indiv.stan'
fit <- stan(modelFile,
            data = data_list,
            iter = 4000,
            warmup = 2000,
            chains = 4,
            cores = 4,
            init =  'random',
            seed = 12345,
            control = list(adapt_delta = 0.9,  max_treedepth = 12)
            )
save(fit, file = paste("fit_subj_",subject,".rdata",sep = ''))

## Plot traces
params <- c('m_min', 'm_plus', 's_min','s_plus')
posterior <- as.array(fit)
np <- nuts_params(fit)
color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior, pars = params, np = np) +
   xlab("Post-warmup iteration")

## summary 
summary(fit, pars = params, , probs = c(0.05, 0.95))$summary

## prediction vs actual
pred <- extract(fit)
choice_pred = sum(pred$y_pred)/(dim(pred$y_pred)[1]*dim(pred$y_pred)[2]) ## predicted proportion of "simultaneous"
choice_pred
mean(data$synchrony) ## actual proportion

##launch_shinystan(fit) ## interface web
