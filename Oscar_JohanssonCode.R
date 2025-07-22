###################### Code used for assignment #############################
# These were all sepearte files which I combined into a single one.
# If something does not work,double check the variables and if the 
# packages are imported.

############### Heavy tail plots ############################
snp <- read.csv('s&p500.csv')
n = 3775
snp$avg <- (snp$Low+snp$High)/2
snp$Date = as.Date.character(snp$Date)
# dates
crash2008 = 'Oct 13 2008'
crash2020 = 'Mar 24 2020'

crashes = snp[which(snp$Date %in% c(crash2008,crash2020)),]

logrtn = log(snp$avg[2:n]/snp$avg[1:(n-1)])

hist(snp$avg[2:n]/snp$avg[1:(n-1)])

which.min(logrtn_stndrd[1:1000])

smu = mean(logrtn)
sigma = sqrt(var(logrtn))
logrtn_stndrd = (logrtn-mu)/sigma
plot(density(logrtn_stndrd,kernel ='gaussian'))

crashes

crashes$x = abs(logrtn_stndrd[as.integer(rownames(crashes))])
crashes$y = 0
crashes$labels = c('2008 crash','2020 crash')

# plot densities
df <- data.frame(
  distribution=factor(rep(c("SNP500", "Normal distribution"), each=n-1)),
  values=c(abs(logrtn_stndrd),abs(rnorm(n-1, mean=0, sd=1)))
)
library(tidyverse)
library(brew)
library(ggplot2)
p<-ggplot(df, aes(x=values, fill=distribution)) +
  geom_density(alpha=0.4,color=NA)+ 
  theme_classic()+
  geom_point(data=crashes,aes(x=x,y=y,label=labels),inherit.aes = FALSE)+
  geom_text(data=crashes,aes(x=x,y=y,label=labels),inherit.aes = FALSE,size=4,vjust=c(-2,-2),hjust=1)+ theme_classic()


p


# QQ plot
qqnorm(logrtn, pch = 1, frame = FALSE)
qqline(logrtn, col = "steelblue", lwd = 2)

install.packages('qqplotr')
gg <- ggplot(data = df, mapping = aes(sample = x)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
gg

# Calculate value at risk
alpha = c(0.1,0.05,0.01,0.001)
var_normal = qnorm(alpha)
cdf = ecdf(logrtn_stndrd)
var_empirical = quantile(cdf,alpha)

########## Define functions ####################################
# Define function
library(EnvStats)
library(spatstat)
library(GoFKernel)
library(miscTools)

## Help Functions
Tfunc <- function(x,a,c){
  M = median(x)
  return( ((x+c)^a-c^a) / ((x+c)^a+(M+c)^a-2*c^a))
}

Tinv <- function(z,a,c,M){
  return( ( (c^a*(2*z-1)-(M+c)^a) / (z-1))^(1/a) -c)
}


l = function(x,a,c){
  N = length(x)
  M = median(x)
  res = N*log(a)
  res = res +N*log((M+c)^a-c^a)
  res = res +(a-1)*sum(log(x+c))
  res = res -2*sum(log((x+c)^a+(M+c)^a-2*c^a))
  return(res)
}

getMLET = function(x){
  ltmp = function(a,c){l(x,a,c)}
  ftmp = function(a){optimize(ltmp, interval=c(0,1000000),a=a,maximum=TRUE)$objective}
  a = optimize(ftmp,interval=c(0.000001,10),maximum=TRUE)$maximum
  c= optimize(ltmp, interval=c(0,1000000),a=a,maximum=TRUE)$maximum
  return(data.frame(c=c,a=a))
}

m <- function(x){
  return(qbeta(0.5*x+0.5,3,3)) #(15/16*(1-x^2)^2)
}

mprim <- function(x){
  return(pbeta(0.5*x+0.5,3,3))#(15/16*(1-2*x^2+x^4))
}

Mfunc <- function(x){
  return(1/16*(15*x-5*x^3+x^5+11))
}


DTKE = function(x,a,method = "w"){
  # Set parameters
  M = median(x)
  n = length(x)
  
  # Get maximum parameters
  res = getMLET(x)
  
  # Get transformed data
  Z = Tfunc(x,res$a,res$c)
  
  # Get beta inverse
  Y = qbeta(Z,3,3)
  Y = Y*2-1
  
  # Calculate bandwiths
  if (method=="w"){
    b = (9/7)^(1/3)*n^(-1/3)
  }else if (method=="x"){
    b=(9/35*m(a)*25/mprim(a))^(1/3)*n^(-1/3)
  }
  
  
  # Get denisity
  DTKE = density(Y,bw=b,kernel="epanechnikov",from=-1,to=1)
  
  # get CDF and inverse
  f =inverse(CDF(DTKE,warn=FALSE))
  q=f(a)
  
  # Calculate VaRa
  VaRa = Tinv(pbeta(q,3,3),res$a,res$c,M)
  return(VaRa)
  
}

TKE = function(x,a){
  # Set parameters
  M = median(x)
  n = length(x)
  
  # Get maximum parameters
  res = getMLET(x)
  
  # Get transformed data
  Z = Tfunc(x,res$a,res$c)
  # Calculate variance
  sigma = sd(Z)
  
  # Rule of thumb AWISE bandwith
  b=sigma*3.572*n^(-1/3)#(8/3)^(1/3)*sigma^(5/3)*n^(-5/3)
  
  # Get denisity
  DTKE = density(Z,bw=b,kernel="epanechnikov",from=0,to=1)
  
  # get CDF and inverse
  f =inverse(CDF(DTKE,warn=FALSE))
  q=f(a)
  
  # Calculate VaRa
  VaRa = Tinv(q,res$a,res$c,M)
  return(VaRa)
  
}

CKE = function(x,a,method = "w"){
  # Set parameters
  sigma = sqrt(var(x))
  mu = mean(x)
  n = length(x)
  x_alpha = qnorm(a)
  
  # Calculate bandwiths
  if (method=="w"){
    b = sigma^(5/3)*(8/3)^(1/3)*n^(-5/3)
  }else if (method=="x"){
    b=(dnorm(x_alpha,mu,sigma)/(ddnorm(x_alpha,mu,sigma)^2)*2/5)^(1/3)*n^(-1/3)
  }
  
  
  # Get denisity
  DTKE = density(x,bw=b,kernel="epanechnikov",from=0)
  
  # get CDF and inverse
  f =inverse(CDF(DTKE,warn=FALSE))
  q=f(a)
  
  # Calculate VaRa
  VaRa = q
  return(VaRa)
  
}

#################### Plot distributions ###################################
n=1000
x= ((0:n))*2/n-1
data = data.frame(x=c(x,x,x),y=c(m(x),M(x),M_false(x)),labels=c(rep('m(x)',n+1),rep('M(x)',n+1),rep('Incorrect M(x)',n+1)))

ggplot(data, aes(x,y,color =labels))+
  geom_line(size=1.5)+
  theme_classic()+ 
  scale_color_manual(values = c("darkred","lightblue", "Royalblue"))

#### Champernowne
T <- function(x,d,c){
  M = median(x)
  return( ((x+c)^d-c^d) / ((x+c)^d+(M+c)^d-2*c^d))
}

T(x,1,5)


x = 1:1000/10
n=5
d = c(0.2,-3,5,10,1)
c = c(0.1,0.5,1,2,5)

chmpData = data.frame(x=x)
for (i in 1:n){
  parameters=sprintf("delta=%.1f,c=%.1f",d[i],c[i])
  chmpData[parameters]<-T(x,d[i],c[i])
}

head(chmpData)
chmpData = chmpData %>% gather(key = "Parameters", value = "y",-x)

ggplot(chmpData, aes(x = x, y = y)) + 
  geom_line(size=0.9,aes(color = Parameters, linetype = Parameters))+
  theme_classic()


################## Evaluated Results as in the report #################
# parameters
iters = 100
n=5000
alpha = 0.99

results = data.frame("DTKE_w"=numeric(),"DTKE_x"=numeric(),"TKE_w"=numeric(),"CKE_w"=numeric(),"CKE_x"=numeric())

for (i in 1:iters){
  
  # Generate data
  x = 0.7*rlnorm(n, 1, 1)+0.3*rpareto(n,1,1)
  cke_w = CKE(x,a=alpha,method="w")
  cke_x = CKE(x,a=alpha,method="x")
  dtke_w = DTKE(x,a=alpha,method="w")
  dtke_x = DTKE(x,a=alpha,method="x")
  tke = TKE(x,a=alpha)
  
  # Add to dataframe
  results[i,] = c(dtke_w,dtke_x,tke,cke_w,cke_x)
}


results %>% gather(key = "Method", value = "esstimate") %>% 
  group_by(Method) %>%
  summarise(mean = mean(esstimate), std = sd(esstimate)) 


##################### Explore Car Claims data ###############################
library(tidyverse)
library(brew)
library(ggplot2)
library(ggplot2)
library(viridis)

# Load data
insurance_data <- read.csv('car_insurance_claims.csv')

# Parameters
alpha = 0.99

# Set variables
claims = insurance_data$CLM_AMT
n = length(claims)
mu = mean(claims)
sigma = sd(claims)

# plot densities
df <- data.frame( Claims = c(claims,rnorm(n,mu,sigma)), 
                  Distribution=c(rep('Car insurance claims',n),rep('lognormal approximation',n)))
ggplot(df, aes(x=Claims,fill = Distribution)) +
  geom_density(alpha=0.6)+
  theme_classic()+
  scale_x_continuous(trans='log10') 



# QQ plot
qqnorm(log(claims), pch = 1, frame = FALSE)
qqline(log(claims), col = "steelblue", lwd = 2)

library(qqplotr)
gg <- ggplot(data = df, mapping = aes(sample = x)) +
  stat_qq_band() +
  stat_qq_line() +
  stat_qq_point() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
gg

# Calculate resulting transformation distributions
res = getMLET(claims)

# Get transformed data
Z = Tfunc(claims,res$a,res$c)

# Plot the fit
empCdf = ecdf(claims)
df = data.frame(Claims = claims, Z=Z, Ecdf=empCdf(claims))
ggplot(df,aes(x=claims))+
  geom_line(aes(y=Z),size=1,color='red')+
  geom_point(aes(y=Ecdf))+
  theme_classic()+
  labs(x = "Claims", y = "Z",title="Fit of Champernowne")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(trans='log10') 


# Get beta inverse
Y = qbeta(Z,3,3)
Y = Y*2-1

# Calculate properties
a=alpha
bw = (9/7)^(1/3)*n^(-1/3)
bx = (9/35*m(a)*25/mprim(a))^(1/3)*n^(-1/3)

# Plot distribution
data1 = data.frame(x=Z)

# Plot distributions
library(ggplot2)
library(viridis)
ggplot(data1) +
  geom_density(alpha=0.5,aes(x=x,fill="bw"),bw=bw,kernel="epanechnikov")+
  geom_density(alpha=0.5,aes(x=x,fill="bx"),bw=bx,kernel="epanechnikov")+
  theme_classic()+
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+
  xlab('Z')+
  ggtitle(sprintf('Density of Z=T(X), a =%0.3f',a))+
  theme(plot.title = element_text(hjust = 0.5))

# Plot distributions'
data2 = data.frame(x=Y)
library(ggplot2)
library(viridis)
ggplot(data2) +
  geom_density(alpha=0.5,aes(x=x,fill="bw"),bw=bw,kernel="epanechnikov")+
  geom_density(alpha=0.5,aes(x=x,fill="bx"),bw=bx,kernel="epanechnikov")+
  theme_classic()+
  scale_color_brewer(palette="Paired")+
  scale_fill_brewer(palette="Paired")+
  xlab('Y')+
  ggtitle(sprintf('Density of Y=M^-1(Z), a =%0.3f',a))+
  theme(plot.title = element_text(hjust = 0.5))


# Evaluate performance
DTKEs = density(Y,bw=bx,kernel="epanechnikov",from=-1,to=1)

# get CDF and inverse
f =inverse(CDF(DTKEs,warn=FALSE))
q=f(a)


# Champernowne
DTKEs = density(Z,bw=bw,kernel="epanechnikov",from=-1,to=1)

# get CDF and inverse
f =inverse(CDF(DTKEs,warn=FALSE))
q=f(a)

################ Evaluate Car claims ########################
# Simulate results
library(EnvStats)
library(dplyr)
library(tidyr)

# parameters
alpha = 0.99

# Load data
insurance_data <- read.csv('car_insurance_claims.csv')

# Set variables
claims = insurance_data$CLM_AMT
# Split
set.seed(1)
n = 5
splits = split(claims, 1:length(claims)%%n)

results = data.frame("DTKE_w"=numeric(),"DTKE_x"=numeric(),"TKE_w"=numeric(),"CKE_w"=numeric(),"CKE_x"=numeric())

for (i in 1:n){
  
  # Generate data
  x = unlist(splits[i])
  cke_w = CKE(x,a=alpha,method="w")
  cke_x = CKE(x,a=alpha,method="x")
  dtke_w = DTKE(x,a=alpha,method="w")
  dtke_x = DTKE(x,a=alpha,method="x")
  tke = TKE(x,a=alpha)
  
  # Add to dataframe
  results[i,] = c(dtke_w,dtke_x,tke,cke_w,cke_x)
}


results %>% gather(key = "Method", value = "esstimate") %>% 
  group_by(Method) %>%
  summarise(mean = mean(esstimate), std = sd(esstimate)) 
