## ---- spaghetti --------

require(rms)
require(data.table)
options(prType='html')    # for model print, summary, anova
getHdata(cdystonia)
setDT(cdystonia)          # convert to data.table

# Construct unique subject ID
cdystonia[, uid := factor(paste(site, id))]

# Tabulate patterns of subjects' time points
cdystonia[, table(tapply(week, uid,
             function(w) paste(sort(unique(w)), collapse=' ')))]

# Plot raw data, superposing subjects
xl <- xlab('Week'); yl <- ylab('TWSTRS-total score')
ggplot(cdystonia, aes(x=week, y=twstrs, color=factor(id))) +
       geom_line() + xl + yl + facet_grid(treat ~ site) +
       guides(color=FALSE)

## ---- quartiles --------

# Show quartiles
cdys <- cdystonia[, j=as.list(quantile(twstrs, (1 : 3)/4)),
                  by = list(treat, week)]
cdys <- upData(cdys, rename=c('25%'='Q1', '50%'='Q2', '75%'='Q3'), print=FALSE)
ggplot(cdys, aes(x=week, y=Q2)) + xl + yl + ylim(0, 70) +
  geom_line() + facet_wrap(~ treat, nrow=2) +
  geom_ribbon(aes(ymin=Q1, ymax=Q3), alpha=0.2)
  
## ---- bootcls --------

# Show means with bootstrap nonparametric CLs
cdys <-  cdystonia[, j=as.list(smean.cl.boot(twstrs)),
                   by = list(treat, week)]
ggplot(cdys, aes(x=week, y=Mean)) + xl + yl + ylim(0, 70) +
  geom_line() + facet_wrap(~ treat, nrow=2) +
  geom_ribbon(aes(x=week, ymin=Lower, ymax=Upper), alpha=0.2)


## ---- e --------

baseline <- cdystonia[week == 0]
baseline[, week := NULL]
setnames(baseline, 'twstrs', 'twstrs0')
followup <- cdystonia[week > 0, .(uid, week, twstrs)]
setkey(baseline, uid)
setkey(followup, uid, week)
both     <- Merge(baseline, followup, id = ~ uid)
# Remove person with no follow-up record
both     <- both[! is.na(week)]
dd       <- datadist(both)
options(datadist='dd')

## ---- k --------

require(nlme)
cp <- list(corCAR1,corExp,corCompSymm,corLin,corGaus,corSpher)
z  <- vector('list',length(cp))
for(k in 1:length(cp)) {
  z[[k]] <- gls(twstrs ~ treat * rcs(week, 3) +
                rcs(twstrs0, 3) + rcs(age, 4) * sex, data=both,
                correlation=cp[[k]](form = ~week | uid))
}
anova(z[[1]],z[[2]],z[[3]],z[[4]],z[[5]],z[[6]])


## ---- l --------

a <- Gls(twstrs ~ treat * rcs(week, 3) + rcs(twstrs0, 3) +
         rcs(age, 4) * sex, data=both,
         correlation=corCAR1(form=~week | uid))
a

## ---- fig-long-variogram --------

v <- Variogram(a, form=~ week | uid)
plot(v)


## ---- fig-long-resid --------

both$resid <- r <- resid(a); both$fitted <- fitted(a)
yl <- ylab('Residuals')
p1 <- ggplot(both, aes(x=fitted, y=resid)) + geom_point() +
      facet_grid(~ treat) + yl
p2 <- ggplot(both, aes(x=twstrs0, y=resid)) + geom_point()+yl
p3 <- ggplot(both, aes(x=week, y=resid)) + yl + ylim(-20,20) +
      stat_summary(fun.data="mean_sdl", geom='smooth')
p4 <- ggplot(both, aes(sample=resid)) + stat_qq() +
      geom_abline(intercept=mean(r), slope=sd(r)) + yl
gridExtra::grid.arrange(p1, p2, p3, p4, ncol=2)


## ---- m --------

anova(a)

## ---- mnote --------

# above was latex(anova(a), file='', label='longit-anova')


## ---- fig-long-anova --------

plot(anova(a))


## ---- fig-long-pleffects --------

ylm <- ylim(25, 60)
p1 <- ggplot(Predict(a, week, treat, conf.int=FALSE),
             adj.subtitle=FALSE, legend.position='top') + ylm
p2 <- ggplot(Predict(a, twstrs0), adj.subtitle=FALSE) + ylm
p3 <- ggplot(Predict(a, age, sex), adj.subtitle=FALSE,
             legend.position='top') + ylm
gridExtra::grid.arrange(p1, p2, p3, ncol=2)


## ---- o --------

summary(a)  # Shows for week 8


## ---- p --------

# To get results for week 8 for a different reference group
# for treatment, use e.g. summary(a, week=4, treat='Placebo')

# Compare low dose with placebo, separately at each time
k1 <- contrast(a, list(week=c(2,4,8,12,16), treat='5000U'),
                  list(week=c(2,4,8,12,16), treat='Placebo'))
options(width=80)
print(k1, digits=3)


## ---- q --------

# Compare high dose with placebo
k2 <- contrast(a, list(week=c(2,4,8,12,16), treat='10000U'),
                  list(week=c(2,4,8,12,16), treat='Placebo'))
print(k2, digits=3)


## ---- fig-long-contrasts --------

k1 <- as.data.frame(k1[c('week', 'Contrast', 'Lower', 'Upper')])
p1 <- ggplot(k1, aes(x=week, y=Contrast)) + geom_point() +
      geom_line() + ylab('Low Dose - Placebo') +
      geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0)
k2 <- as.data.frame(k2[c('week', 'Contrast', 'Lower', 'Upper')])
p2 <- ggplot(k2, aes(x=week, y=Contrast)) + geom_point() +
      geom_line() + ylab('High Dose - Placebo') +
      geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0)
gridExtra::grid.arrange(p1, p2, ncol=2)




## ---- nomogram --------

n <- nomogram(a, age=c(seq(20, 80, by=10), 85))
plot(n, cex.axis=.55, cex.var=.8, lmgp=.25)


## ---- bayesfit --------

require(rmsb)
cmdstanr::set_cmdstan_path(cmdstan.loc)
# cmdstan.loc is defined in ~/.Rprofile
options(mc.cores=parallel::detectCores() - 1, rmsb.backend='cmdstan')
bpo <- blrm(twstrs ~ treat * rcs(week, 3) + rcs(twstrs0, 3) +
            rcs(age, 4) * sex + cluster(uid), data=both, file='bpo.rds')
# file= means that after the first time the model is run, it will not
# be re-run unless the data, fitting options, or underlying Stan code change
stanDx(bpo)

## ---- bayesfit2 --------

print(bpo, intercepts=TRUE)
a <- anova(bpo)
a
plot(a)

## ---- bayesfit3 --------

wks <- c(2,4,8,12,16)
k <- contrast(bpo, list(week=wks, treat='10000U'),
                   list(week=wks, treat='Placebo'),
              cnames=paste('Week', wks))
k
plot(k)

## ---- bayesfit4 --------

k <- as.data.frame(k[c('week', 'Contrast', 'Lower', 'Upper')])
ggplot(k, aes(x=week, y=Contrast)) + geom_point() +
  geom_line() + ylab('High Dose - Placebo') +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0)

## ---- bayesfit5 --------

M <- Mean(bpo)   # create R function that computes mean Y from X*beta
k <- contrast(bpo, list(week=wks, treat='10000U'),
                   list(week=wks, treat='Placebo'),
              fun=M, cnames=paste('Week', wks))
plot(k, which='diff') + theme(legend.position='bottom')

## ---- bayesfit6 --------

f <- function(x) {
  hpd <- HPDint(x, prob=0.95)   # is in rmsb
  r <- c(mean(x), median(x), hpd)
  names(r) <- c('Mean', 'Median', 'Lower', 'Upper')
  r
}
w    <- as.data.frame(t(apply(k$esta - k$estb, 2, f)))
week <- as.numeric(sub('Week ', '', rownames(w)))
ggplot(w, aes(x=week, y=Mean)) + geom_point() +
  geom_line() + ylab('High Dose - Placebo') +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0) +
  scale_y_continuous(breaks=c(-8, -4, 0, 4))

## ---- mark1 --------

# Create a new variable to hold previous value of Y for the subject
# For week 2, previous value is the baseline value
setDT(both, key=c('uid', 'week'))
both[, ptwstrs := shift(twstrs), by=uid]
both[week == 2, ptwstrs := twstrs0]
dd <- datadist(both)
bmark <- blrm(twstrs ~  treat * rcs(week, 3) + rcs(ptwstrs, 4) +
                        rcs(age, 4) * sex,
              data=both, file='bmark.rds')
# When adding partial PO terms for week and ptwstrs, z=-1.8, 5.04
stanDx(bmark)
stanDxplot(bmark)


## ---- mark2 --------

bmark
a <- anova(bmark)
a
plot(a)

## ---- mark3 --------

bmarkre <- blrm(twstrs ~  treat * rcs(week, 3) + rcs(ptwstrs, 4) +
                          rcs(age, 4) * sex + cluster(uid),
                data=both, file='bmarkre.rds')
stanDx(bmarkre)

## ---- mark4 --------

bmarkre

## ---- mark4r --------

plot(sqrt(diag(vcov(bmark))), sqrt(diag(vcov(bmarkre))),
     xlab='Posterior SDs in Conditional Independence Markov Model',
     ylab='Posterior SDs in Random Effects Markov Model')
abline(a=0, b=1, col=gray(0.85))

## ---- mark5 --------

ggplot(Predict(bmark))
ggplot(Predict(bmark, week, treat))

## ---- mark5b --------

k <- contrast(bmark, list(week=wks, treat='10000U'),
                     list(week=wks, treat='Placebo'),
              cnames=paste('Week', wks))
k
plot(k)
k <- as.data.frame(k[c('week', 'Contrast', 'Lower', 'Upper')])
ggplot(k, aes(x=week, y=Contrast)) + geom_point() +
  geom_line() + ylab('High Dose - Placebo') +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=0)

## ---- mark6 --------

ex <- ExProb(bmark)
ex40 <- function(lp, ...) ex(lp, y=40, ...)
ggplot(Predict(bmark, week, treat, ptwstrs=40, fun=ex40))

## ---- mark6b --------

ggplot(Predict(bmark, week, treat, ptwstrs=40, fun=Mean(bmark)))

## ---- mark6c --------

# Get median/mode for covariates including ptwstrs (TWSTRS in previous visit)
d <- gendata(bmark)
d
d$week <- 12
p <- predict(bmark, d, type='fitted.ind')   # defaults to posterior means
yvals <- as.numeric(sub('twstrs=', '', p$y))
lo <- p$Lower / sum(p$Lower)
hi <- p$Upper / sum(p$Upper)
plot(yvals, p$Mean, type='l', xlab='TWSTRS', ylab='',
     ylim=range(c(lo, hi)))
lines(yvals, lo, col=gray(0.8))
lines(yvals, hi, col=gray(0.8))

## ---- mark6d --------

p <- predict(bmark, d, type='fitted.ind', posterior.summary='all')
cols <- adjustcolor(1 : 10, 0.7)
for(i in 1 : 5) {
  if(i == 1) plot(yvals, p[i, 1, ], type='l', col=cols[1], xlab='TWSTRS', ylab='')
  else lines(yvals, p[i, 1, ], col=cols[i])
}

## ---- mark7 --------

# Baseline twstrs to 42 in d
# For each dose, get all the posterior draws for all state occupancy
# probabilities for all visit
ylev <- sort(unique(both$twstrs))
tlev <- c('Placebo', '10000U')
R <- list()
for(trt in tlev) {   # separately by treatment
  d$treat <- trt
  u <- soprobMarkovOrdm(bmark, d, wks, ylev,
                        tvarname='week', pvarname='ptwstrs')
  R[[trt]] <- u
}
dim(R[[1]])    # posterior draws x times x distinct twstrs values

# For each posterior draw, treatment, and week compute the mean TWSTRS
# Then compute posterior mean of means, and HPD interval
Rmean <- Rmeans <- list()
for(trt in tlev) {
  r <- R[[trt]]
  # Mean Y at each week and posterior draw (mean from a discrete distribution)
  m <- apply(r, 1:2, function(x) sum(ylev * x))
  Rmeans[[trt]] <- m
  # Posterior mean and median and HPD interval over draws
  u <- apply(m, 2, f)   # f defined above
  u <- rbind(week=as.numeric(colnames(u)), u)
  Rmean[[trt]] <- u
}
r <- lapply(Rmean, function(x) as.data.frame(t(x)))
for(trt in tlev) r[[trt]]$treat <- trt
r <- do.call(rbind, r)
ggplot(r, aes(x=week, y=Mean, color=treat)) + geom_line() +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2, linetype=0)

## ---- mark8 --------

Dif <- Rmeans$`10000U` - Rmeans$Placebo
dif <- as.data.frame(t(apply(Dif, 2, f)))
dif$week <- as.numeric(rownames(dif))
ggplot(dif, aes(x=week, y=Mean)) + geom_line() +
  geom_ribbon(aes(ymin=Lower, ymax=Upper), alpha=0.2, linetype=0) +
  ylab('High Dose - Placebo TWSTRS')

## ---- mark9 --------

p <- R$`10000U`[, '12', ]   # 4000 x 62
pmean <- apply(p, 2, mean)
yvals <- as.numeric(names(pmean))
plot(yvals, pmean, type='l', xlab='TWSTRS', ylab='')

