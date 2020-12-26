library(ecoCopula)
library(mvabund)

data(spider)
spiddat = mvabund(spider$abund)

pa = spiddat
pa[pa>0] = 1

abund = spider$abund
abund[,1:3] = (abund[,1:3]>0)*1
myfamily = c(
  rep(c("binomial"), 3),
  rep(c("negative.binomial"), (ncol(abund)-3))
)
abund1 = spider$abund

X = data.frame(spider$x)
X$Treatment2 = rep(c("A","B"),each=14)
X$Treatment4 = rep(c("A","B","C","D"),each=7)
X$month = rep(c("Jan","Feb","Mar","Apr"),each=7)

Xnew = X
Xnew$Treatment2[6:7] = c("B","B")
Xnew_sub = Xnew[1:10,]

Xvec = data.frame(Treatment4 = X$Treatment4[1:15])

fit1.lm = manylm(spiddat~bare.sand, data=X)
fit1.cord = cord(fit1.lm)

nb0.glm = manyglm(spiddat~1, family="negative.binomial")
nb0.cord = cord(nb0.glm)

nb0.ssdm = stackedsdm(abund1, ~1, data=X, family="negative.binomial", ncores=2)
nb0.ssdm.cord = cord(nb0.ssdm)

nb1.ssdm = stackedsdm(abund1, ~soil.dry, data=X, family="negative.binomial", ncores=2)
nb1.ssdm.cord = cord(nb1.ssdm)

ssdm_multi = stackedsdm(abund, formula_X = ~ bare.sand, data = X, family = myfamily, ncores = 2)
ssdm_multi.cord = cord(ssdm_multi)

nb1.glm = manyglm(spiddat~bare.sand, family="negative.binomial", data=X)
nb1.cord = cord(nb1.glm)

bin0.glm = manyglm(pa~1, family="binomial")
bin0.cord = cord(bin0.glm)

poi1.glm = manyglm(spiddat~bare.sand, family="poisson", data=X)
poi1.cord = cord(poi1.glm)

nb2.glm = manyglm(spiddat~soil.dry+bare.sand, family="negative.binomial", data=X)
nb2.cord = cord(nb2.glm)

nb_mix.glm = manyglm(spiddat~month+bare.sand, family="negative.binomial", data=X)
nb_mix.cord = cord(nb_mix.glm)

nb_fac2.glm = manyglm(spiddat~Treatment2, family="negative.binomial", data=X)
nb_fac2.cord = cord(nb_fac2.glm)

nb_fac4.glm = manyglm(spiddat~Treatment4, family="negative.binomial", data=X)
nb_fac4.cord = cord(nb_fac4.glm)

nb_fac4.ssdm = stackedsdm(abund1, ~Treatment4, data=X, family="negative.binomial", ncores=2)
nb_fac4.ssdm.cord = cord(nb_fac4.ssdm)

nb_mth.glm = manyglm(spiddat~month+Treatment2, family="negative.binomial", data=X)
nb_mth.cord = cord(nb_mth.glm)
