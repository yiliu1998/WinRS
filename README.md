# WinRS
`WinRS` is an R package that implements the nonparametric maximum likelihood estimation (NPMLE) of the win ratio (WR) for two hierarchical endpoints in a clinical study with a finite time horizon. We focus on the setting where the first endpoint is a terminal time-to-event outcome subject to right-censoring, and the second is a non-survival endpoint—measured only once at the end of study (e.g., a quality-of-life or response outcome)—subject to missingness. The NPMLE is based on our proposed method called "S score" in the paper, also referred to as the "S score estimator." 

## Installation
To install the latest version of the R package from GitHub, please run following commands in R:

```r
if (!require("devtools"))
install.packages("devtools")
devtools::install_github("yiliu1998/WinRS")
```

## Usage
Here is a toy example for using the package. 
```r
# function for generating the toy example data
dat_gen_cov <- function(n=10000,
                        hazard_T=0.02,hazard_C=0.01,h=90,censoring=T,
                        mu=10,sigma=5,tau=50,beta=c(-1,0.5,0.8,0.3)) {

  X1 <- rnorm(n)
  X2 <- rbinom(n,1,0.5)
  X3 <- runif(n,-1,1)

  T1 <- rexp(n,rate=hazard_T) %>% round()
  if(censoring) {
    C1 <- rexp(n,rate=hazard_C) %>% round()
  } else {
    C1 <- h+1
  }

  Y1 <- pmin(T1,h)+I(T1>h)
  Y1.tilde <- pmin(Y1,C1,h)+I(Y1>h & C1>h)
  Delta1 <- as.numeric(Y1<=C1 & Y1<=h)+as.numeric(Y1>h & C1>h)

  Y2 <- (pmin(pmax(rnorm(n,mean=mu,sd=sigma),0),tau)*I(Y1>h)) %>% round()
  logit_p2 <- beta[1]+beta[2]*X1+beta[3]*X2+beta[4]*X3
  p2 <- 1 / (1 + exp(-logit_p2)) * I(Y1>h)
  R2 <- runif(n)>=p2
  Y2[R2==0] <- NA

  S <- Y1.tilde+as.numeric(Y1.tilde>h)*ifelse(R2==0, 0, Y2)
  Delta_S <- ifelse((Delta1==1 & Y1.tilde<=h) | (R2==1 & Y1.tilde>h), 1, 0)

  data.frame(X1=X1, X2=X2, X3=X3,
             Y1=Y1.tilde, Delta1=Delta1, Y2=Y2, R2=as.numeric(R2),
             S=S, Delta_S=Delta_S)
}

# generate data of two treatment groups
df_a <- dat_gen_cov(n=10000,
                      hazard_T=0.02, hazard_C=0.01, h=90, censoring=T,
                      mu=10, sigma=5, tau=50, beta=c(-1,0.5,0.8,0.3))
df_b <- dat_gen_cov(n=10000,
                      hazard_T=0.03, hazard_C=0.01, h=90, censoring=T,
                      mu=10, sigma=10, tau=50, beta=c(-1,0.5,0.8,0.3))
h=90
covar=c("X1", "X2", "X3")
WinRS(df_a, df_b,
      cov.adjust=TRUE, covar.name=covar,
      t.horizon=h, out1.name="Y1", delta1.name="Delta1",
      out2.name="Y2",
      boot=TRUE, n_boot=100)
```

The following is the output of the above program:
```r
$WR.est
[1] 1.486372

$IF.sd
[1] 0.03129546

$IF.lwr
[1] 1.425033

$IF.upr
[1] 1.547711

$boot.sd
[1] 0.02524311

$boot.Wald.lwr
[1] 1.436895

$boot.Wald.upr
[1] 1.535848

$boot.qt.lwr
    2.5% 
1.437287 

$boot.qt.upr
   97.5% 
1.537708 
```

## Reference
Please cite the following paper: 

TBA... 

## Contact
Please reach out to Yi Liu (yi.liu.biostat@gmail.com) if you have any questions for the package. 
