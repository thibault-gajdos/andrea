rm(list=ls(all=TRUE))  ## efface les données
source('~/thib/projects/tools/R_lib.r')
setwd('~/thib/projects/andrea_stan/')

## Caution! times (in seconds) are multiplied by 10 for computational issues. Thus, output must be multiplied by 100 to get results in ms

## * Prepare data 
data <- read_csv('ares_expt1_bhv.csv') %>%
    mutate(SOAactual =  SOAactual*10) %>%
    group_by(subject) %>%
    mutate(trial =  row_number())%>%
    ungroup()

## build data list
N <- length(unique(data$subject))
T <- length(unique(data$trial))
subjs <- unique(data$subject)

## compute trials by subject
d <- data %>%
    group_by(subject) %>%
    summarise(t_subjs = n(), std_dev = sd(SOAactual)) 
t_subjs <- d$t_subjs
std_dev <- d$std_dev
## Initialize data arrays
delay    <- array(0, c(N, T))
choice <- array(0, c(N, T))

## Write from raw_data to the data arrays
for (i in 1:N) {
    t <- t_subjs[i]
    data_subj <- data %>% filter(subject == subjs[i])
    delay[i, 1:t]    <- data_subj$SOAactual
    choice[i, 1:t] <- data_subj$synchrony
}

## Wrap into a list for Stan
data_list <- list(
    N        = N,
    T        = T,
    Tsubj    = t_subjs,
    std_dev = std_dev,
    choice   = choice,
    delay  = delay
    )

## * fit
modelFile <- '~/thib/projects/andrea_stan/asym_hierarchical.stan'
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
save(fit, file = 'fit_hierarchical.rdata')

## Plot traces
params <- c('mu','sigma')
posterior <- as.array(fit)
np <- nuts_params(fit)
color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior, pars = params, np = np) +
   xlab("Post-warmup iteration")


## summary
## group level parameters
params <- c('mu','sigma')
summary(fit, pars = params, , probs = c(0.05, 0.95))$summary


## |         |       mean|   se_mean|        sd|         5%|        95%|    n_eff|      Rhat|
## |:--------|----------:|---------:|---------:|----------:|----------:|--------:|---------:|
## |mu[1]    | -0.7954720| 0.0032774| 0.1264954| -1.0027020| -0.5911700| 1489.627| 1.0016429|
## |mu[2]    |  4.0141046| 0.0066389| 0.2775774|  3.5697909|  4.4831634| 1748.113| 1.0038406|
## |mu[3]    | -2.0756775| 0.0064112| 0.4991847| -2.9906878| -1.3737949| 6062.452| 1.0002914|
## |mu[4]    | -1.6467289| 0.0063777| 0.5596160| -2.6598692| -0.8297820| 7699.336| 0.9996907|
## |sigma[1] |  0.4888036| 0.0022790| 0.1060312|  0.3428116|  0.6827371| 2164.588| 1.0006047|
## |sigma[2] |  1.0098712| 0.0044889| 0.2169379|  0.7075190|  1.4098364| 2335.577| 1.0069989|
## |sigma[3] |  0.3010226| 0.0030205| 0.2530945|  0.0202779|  0.7998223| 7021.114| 0.9997656|
## |sigma[4] |  0.3643663| 0.0034000| 0.3053866|  0.0229669|  0.9613414| 8067.347| 0.9999273|


## individual level parameters
params <- c('m_plus','m_min')
summary(fit, pars = params, , probs = c(0.05, 0.95))$summary

## |           |       mean|   se_mean|        sd|         5%|        95%|     n_eff|      Rhat|
## |:----------|----------:|---------:|---------:|----------:|----------:|---------:|---------:|
## |m_plus[1]  |  2.5160393| 0.0021285| 0.2099433|  2.1776231|  2.8676511|  9728.544| 0.9999816|
## |m_plus[2]  |  3.1017782| 0.0020656| 0.2086368|  2.7640930|  3.4473851| 10202.275| 0.9999741|
## |m_plus[3]  |  4.2472057| 0.0026696| 0.3027430|  3.7729263|  4.7636707| 12860.168| 0.9995578|
## |m_plus[4]  |  4.0292275| 0.0028821| 0.3053485|  3.5524271|  4.5476276| 11224.710| 1.0001119|
## |m_plus[5]  |  5.2458809| 0.0046805| 0.4672765|  4.5414151|  6.0539945|  9967.131| 0.9999253|
## |m_plus[6]  |  3.9973759| 0.0027748| 0.2972435|  3.5295508|  4.5046827| 11475.124| 1.0000455|
## |m_plus[7]  |  4.7972113| 0.0050956| 0.4653051|  4.1374920|  5.6093610|  8338.539| 0.9999287|
## |m_plus[8]  |  4.4048777| 0.0028106| 0.3221427|  3.8965511|  4.9554764| 13137.118| 1.0000922|
## |m_plus[9]  |  3.7882137| 0.0021515| 0.2464686|  3.4003122|  4.2066568| 13123.190| 0.9996397|
## |m_plus[10] |  2.4687426| 0.0020228| 0.1977085|  2.1540058|  2.8057838|  9552.682| 1.0003398|
## |m_plus[11] |  4.8393506| 0.0037118| 0.3912538|  4.2484535|  5.5104272| 11111.045| 0.9998701|
## |m_plus[12] |  2.7220130| 0.0019258| 0.1968789|  2.4056083|  3.0558816| 10450.915| 0.9998764|
## |m_plus[13] |  4.3702085| 0.0036571| 0.3674295|  3.8236312|  5.0122949| 10094.390| 1.0002919|
## |m_plus[14] |  4.5764996| 0.0040617| 0.4170841|  3.9396663|  5.3087626| 10544.735| 0.9996157|
## |m_plus[15] |  3.9746274| 0.0024630| 0.2788002|  3.5355893|  4.4437531| 12812.925| 0.9998850|
## |m_plus[16] |  5.2988853| 0.0083643| 0.6598206|  4.4100010|  6.5234815|  6222.881| 1.0013264|
## |m_min[1]   | -0.6032446| 0.0012596| 0.1405812| -0.8361374| -0.3745196| 12456.473| 0.9997941|
## |m_min[2]   | -1.3261882| 0.0012265| 0.1346075| -1.5468415| -1.1092603| 12044.284| 0.9997369|
## |m_min[3]   | -1.5188876| 0.0012857| 0.1352777| -1.7447693| -1.3000884| 11069.910| 0.9998366|
## |m_min[4]   | -0.5764175| 0.0011427| 0.1258728| -0.7878682| -0.3706188| 12133.802| 0.9996895|
## |m_min[5]   | -1.4462042| 0.0013394| 0.1422899| -1.6817662| -1.2172012| 11285.270| 1.0000442|
## |m_min[6]   | -0.2987280| 0.0011071| 0.1254947| -0.5076135| -0.0930675| 12848.267| 0.9996731|
## |m_min[7]   | -0.7462670| 0.0011515| 0.1244078| -0.9530753| -0.5409995| 11672.527| 0.9998116|
## |m_min[8]   | -0.7060730| 0.0011176| 0.1309795| -0.9191227| -0.4902540| 13734.403| 0.9996274|
## |m_min[9]   | -0.1780393| 0.0012483| 0.1288644| -0.3926115|  0.0324752| 10657.553| 0.9997689|
## |m_min[10]  | -0.7697682| 0.0011813| 0.1383698| -0.9977515| -0.5448672| 13721.172| 0.9998271|
## |m_min[11]  | -0.9823562| 0.0011177| 0.1236892| -1.1865277| -0.7794518| 12247.034| 0.9997366|
## |m_min[12]  | -0.7802530| 0.0010118| 0.1168230| -0.9746885| -0.5877001| 13330.880| 1.0001560|
## |m_min[13]  | -0.2436710| 0.0011929| 0.1210090| -0.4449507| -0.0469946| 10290.281| 1.0000681|
## |m_min[14]  | -1.4309309| 0.0012300| 0.1265212| -1.6415198| -1.2220205| 10581.479| 0.9996677|
## |m_min[15]  | -0.4483699| 0.0010710| 0.1203661| -0.6469546| -0.2534818| 12629.640| 1.0002939|
## |m_min[16]  | -0.6155175| 0.0010617| 0.1219662| -0.8144302| -0.4162011| 13196.495| 0.9998814|



## prediction vs actual
## produce a dataframe with actual and simulated proportion of synchrony judgement by subject
d <- data %>%
    group_by(subject) %>%
    summarise(t_subjs = n(), synchrony = mean(synchrony))  
t <- d$t_subjs
synchrony = d$synchrony
pred = data.frame(subject = numeric(),synchrony = numeric(),sim_synchrony = numeric())
for (i in c(1:N)){
    pred[i,1] <- i
    pred[i,2] <- synchrony[i]
    pred[i,3] <- table(extract(fit)$y_pred[,i,1:t[i]])[2]/(t[i]*dim(extract(fit)$y_pred[,i,1:t[i]])[1])
}
save(pred, file = 'pred.rdata')
## | subject| synchrony| sim_synchrony|
## |-------:|---------:|-------------:|
## |       1| 0.4642857|     0.4751101|
## |       2| 0.6398810|     0.6474967|
## |       3| 0.7172619|     0.7162269|
## |       4| 0.5446429|     0.5506682|
## |       5| 0.6964286|     0.6911737|
## |       6| 0.4910714|     0.5036507|
## |       7| 0.5952381|     0.6004944|
## |       8| 0.5416667|     0.5490104|
## |       9| 0.4880952|     0.5060807|
## |      10| 0.4910714|     0.5056302|
## |      11| 0.6517857|     0.6552180|
## |      12| 0.5267857|     0.5451693|
## |      13| 0.5059524|     0.5138597|
## |      14| 0.6726190|     0.6691239|
## |      15| 0.5084337|     0.5168346|
## |      16| 0.5714286|     0.5731544|

## Draws group parameters (m_plus, m_min, s_min, s_plus)
dd <- as.data.frame(extract(fit)$mu) %>%
    rename(m_min = V1, m_plus = V2, s_min = V3, s_plus = V4)
head(dd)


launch_shinystan(fit) ## interface web
