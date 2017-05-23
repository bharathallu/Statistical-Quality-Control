rm(list=ls())
data<- read.csv("C:/ISEN 614/Project/proj.csv",header = F,sep = ",")
summary(data)
data_scale<- scale(data,center = T,scale = T)

### t2 square for the data without compression#####(need for dimensional reduction)
MU <- colMeans(data)
VAR <- var(data)
DIST <- c(rep(0,552))

for (i in 1:552)
  {
  DIST[i] <- as.matrix(data[i,]-MU) %*% solve(as.matrix(VAR)) %*% t(as.matrix(data[i,]-MU))
}

par(mfrow=c(1,1))
plot(DIST,type="b",main = "Hotelling T2: Individual Observations",ylab="T2 statistic",xlab="observations")
abline(h=qchisq(1-.0027,209),col="red",lwd=3,lty=2)

sum(DIST > qchisq(1-.0027,209))

qchisq(1-.0027,209)


data_pca <- princomp(data)
data_pca_scale<- princomp(data_scale)
names(data_pca)


as <- eigen(var(data))
su <- sum(as$values)
prop <- as$values/su
barplot(prop,col="green",ylim=c(0,.7),xlim=c(0,100),ylab = "proportion of variance explained",xlab = "index",main="Pareto Plot for Variance by each PC")


par(mfrow=c(1,1))
screeplot(data_pca,npcs=80,type = "lines",main="Screeplot")
screeplot(data_pca_scale,npcs=80,type="lines")


PCA1<- as.matrix(data)%*%data_pca_scale$loadings[,1]
PCA2<- as.matrix(data)%*%data_pca_scale$loadings[,2]
PCA <- as.matrix(data)%*%data_pca$loadings[,1:3]
par(mfrow=c(1,1))
plot(PCA1,PCA2)

mu<- colMeans(data)
sigma<- var(data_scale)
library(mvtnorm)
X <- rmvnorm(1000, mu, sigma)
library(energy)
mvnorm.etest(X)


dim(mu)
dim(sigma)
head(data_scale)

eig <- eigen(var(data))

p <- length(data)
n <- 509

mdl <- c(rep(0,209))
a_l <- c(rep(0,209))
g_l <- c(rep(0,209))

for (l in 0:208) 
  {
  a_l[l] <- mean(eig$values[1:(p-l)])
  g_l[l] <- exp(mean(log(eig$values[1:(p-l)])))
 
}

for (l in 0:208) {
  mdl[l] <- (n*(p-l)*log(a_l[l]/g_l[l]))
  
}

####uni-variate control charts after pca ####

PCA <- as.matrix(data)%*%data_pca$loadings[,1:3]

plot(PCA[,1],type="b",main="Univariate X bar chart for PCA1",ylab="PCA1")
abline(h= mean(PCA[,1]),col="blue",lwd=3,lty=2)
abline(h= mean(PCA[,1]) + (qnorm(1-.0009/2)*sqrt(var(PCA[,1]))),col="red",lwd=3,lty=2)
abline(h= mean(PCA[,1]) - (qnorm(1-.0009/2)*sqrt(var(PCA[,1]))),col="red",lwd=3,lty=2)

plot(PCA[,2],type="b",main="Univariate X bar chart for PCA2",ylab="PCA2",ylim=c(50,500))
abline(h= mean(PCA[,2]),col="blue",lwd=3,lty=2)
abline(h= mean(PCA[,2]) + (qnorm(1-.0009/2)*sqrt(var(PCA[,2]))),col="red",lwd=3,lty=2)
abline(h= mean(PCA[,2]) - (qnorm(1-.0009/2)*sqrt(var(PCA[,2]))),col="red",lwd=3,lty=2)

plot(PCA[,3],type="b",main="Univariate X bar chart for PCA3",ylab="PCA3",ylim=c(500,800))
abline(h= mean(PCA[,3]),col="blue",lwd=3,lty=2)
abline(h= mean(PCA[,3]) + (qnorm(1-.0009/2)*sqrt(var(PCA[,3]))),col="red",lwd=3,lty=2)
abline(h= mean(PCA[,3]) - (qnorm(1-.0009/2)*sqrt(var(PCA[,3]))),col="red",lwd=3,lty=2)


##### multi variate t2 chart#####

PCA <- as.matrix(data)%*%data_pca$loadings[,1:3]
mu <- colMeans(PCA)
var_pca <- var(PCA)
dist <- c(rep(0,552))

for (i in 1:552){
dist[i] <- t(PCA[i,]-mu) %*% solve(var_pca) %*% (PCA[i,]-mu)
}

par(mfrow=c(1,1))
plot(dist,type="b",main="Hotelling T2 on Reduced Data: 1st Iteration",ylab="T2")
abline(h=qchisq(1-.0027,3),col="red",lwd=3,lty=2)

dist[dist > qchisq(1-.0027,3)]
sum(dist > qchisq(1-.0027,3))

x1 <- PCA[which(dist > qchisq(1-.0027,3)),]




PCA1 <- PCA[dist < qchisq(1-.0027,3),]
mu1 <- colMeans(PCA1)
var_pca1 <- var(PCA1)
dist1 <- c(rep(0,nrow(PCA1)))

for (i in 1:nrow(PCA1)){
  dist1[i] <- t(PCA1[i,]-mu1) %*% solve(var_pca1) %*% (PCA1[i,]-mu1)
}

dist1[dist1 > qchisq(1-.0027,3)]
sum(dist1 > qchisq(1-.0027,3))

par(mfrow=c(1,1))
plot(dist1,type="b",main="Hotelling T2 on Reduced Data: 2nd Iteration",ylab="T2")
abline(h=qchisq(1-.0027,3),col="red",lwd=3,lty=2)

x2 <- PCA1[which(dist1 > qchisq(1-.0027,3)),]

PCA2 <- PCA1[dist1 < qchisq(1-.0027,3),]
  mu2 <- colMeans(PCA2)
var_pca2 <- var(PCA2)
dist2 <- c(rep(0,nrow(PCA2)))

for (i in 1:nrow(PCA2)){
  dist2[i] <- t(PCA2[i,]-mu2) %*% solve(var_pca2) %*% (PCA2[i,]-mu2)
}

dist2[dist2 > qchisq(1-.0027,3)]
sum(dist2 > qchisq(1-.0027,3))
x3 <- PCA2[which(dist2 > qchisq(1-.0027,3)),]

##### Iteration 3 ###########

PCA3 <- PCA2[dist2 < qchisq(1-.0027,3),]
mu3 <- colMeans(PCA3)
var_pca3 <- var(PCA3)
dist3 <- c(rep(0,nrow(PCA3)))

for (i in 1:nrow(PCA3)){
  dist3[i] <- t(PCA3[i,]-mu3) %*% solve(var_pca3) %*% (PCA3[i,]-mu3)
}

dist3[dist3 > qchisq(1-.0027,3)]
sum(dist3 > qchisq(1-.0027,3))
x4 <- PCA3[which(dist3 > qchisq(1-.0027,3)),]

###### Iteration 4 ###########

PCA4 <- PCA3[dist3 < qchisq(1-.0027,3),]
mu4 <- colMeans(PCA4)
var_pca4 <- var(PCA4)
dist4 <- c(rep(0,nrow(PCA4)))

for (i in 1:nrow(PCA4)){
  dist4[i] <- t(PCA4[i,]-mu4) %*% solve(var_pca4) %*% (PCA4[i,]-mu4)
}

dist4[dist4 > qchisq(1-.0027,3)]
sum(dist4 > qchisq(1-.0027,3))
x5 <- PCA4[which(dist4 > qchisq(1-.0027,3)),]

##### Iteration 5 ######

PCA5 <- PCA4[dist4 < qchisq(1-.0027,3),]
mu5 <- colMeans(PCA5)
var_pca5 <- var(PCA5)
dist5 <- c(rep(0,nrow(PCA5)))

for (i in 1:nrow(PCA5)){
  dist5[i] <- t(PCA5[i,]-mu5) %*% solve(var_pca5) %*% (PCA5[i,]-mu5)
}

dist5[dist5 > qchisq(1-.0027,3)]
sum(dist5 > qchisq(1-.0027,3))
x6 <- PCA5[which(dist5 > qchisq(1-.0027,3)),]
###### Iteration 6 #######

PCA6 <- PCA5[dist5 < qchisq(1-.0027,3),]
mu6 <- colMeans(PCA6)
var_pca6 <- var(PCA6)
dist6 <- c(rep(0,nrow(PCA6)))

for (i in 1:nrow(PCA6)){
  dist6[i] <- t(PCA6[i,]-mu6) %*% solve(var_pca6) %*% (PCA6[i,]-mu6)
}

dist6[dist6 > qchisq(1-.0027,3)]
sum(dist6 > qchisq(1-.0027,3))
x7 <- PCA6[which(dist6 > qchisq(1-.0027,3)),]
###### Iteration 7 #########

PCA7 <- PCA6[dist6 < qchisq(1-.0027,3),]
mu7 <- colMeans(PCA7)
var_pca7 <- var(PCA7)
dist7 <- c(rep(0,nrow(PCA7)))

for (i in 1:nrow(PCA7)){
  dist7[i] <- t(PCA7[i,]-mu7) %*% solve(var_pca7) %*% (PCA7[i,]-mu7)
}

dist7[dist7 > qchisq(1-.0027,3)]
sum(dist7 > qchisq(1-.0027,3))
x8 <- PCA7[which(dist7> qchisq(1-.0027,3)),]
####### Iteration 8 ######

PCA8 <- PCA7[dist7 < qchisq(1-.0027,3),]
mu8 <- colMeans(PCA8)
var_pca8 <- var(PCA8)
dist8 <- c(rep(0,nrow(PCA8)))

for (i in 1:nrow(PCA8)){
  dist8[i] <- t(PCA8[i,]-mu8) %*% solve(var_pca8) %*% (PCA8[i,]-mu8)
}

dist8[dist8 > qchisq(1-.0027,3)]
sum(dist8 > qchisq(1-.0027,3))
x9 <- PCA8[which(dist8 > qchisq(1-.0027,3)),]
## Iteration 9 ##
PCA9 <- PCA8[dist8 < qchisq(1-.0027,3),]
mu9 <- colMeans(PCA9)
var_pca9 <- var(PCA9)
dist9 <- c(rep(0,nrow(PCA9)))

for (i in 1:nrow(PCA9)){
  dist9[i] <- t(PCA9[i,]-mu9) %*% solve(var_pca9) %*% (PCA9[i,]-mu9)
}

dist9[dist9 > qchisq(1-.0027,3)]
sum(dist9 > qchisq(1-.0027,3))
PCA9[which(dist9 > qchisq(1-.0027,3)),]

par(mfrow=c(1,1))
plot(dist9,type="b",main="Hotelling T2 on Reduced Data: 10th Iteration",ylab="T2",ylim=c(-2,16))
abline(h=qchisq(1-.0027,3),col="red",lwd=3,lty=2)


########mcusum after t2 ##########

dta_cusum.mqcdt <- mqcd(PCA9)
res_t <- mqcs.mcusum(dta_cusum.mqcdt,h=6.25,k=1.5,plot=T)
summary(res_t)           
res_t$violations
x_10 <- PCA9[res_t$violations,]

PCA12 <- PCA9[-res_t$violations,]
dta_cusum.mqcd1 <- mqcd(PCA12)
res1 <- mqcs.mcusum(dta_cusum.mqcd1,h=6.25,k=1.5,plot=T)
summary(res1)           
res1$violations
x_11 <- PCA12[res1$violations,]

PCA13 <- PCA12[-res1$violations,]
dta_cusum.mqcd2 <- mqcd(PCA13)
res2 <- mqcs.mcusum(dta_cusum.mqcd2,h=6.25,k=1.5,plot=T)
summary(res2)           
res2$violations
x_12 <- PCA13[res2$violations,]

PCA14 <- PCA13[-res2$violations,]
dta_cusum.mqcd3 <- mqcd(PCA14)
res3 <- mqcs.mcusum(dta_cusum.mqcd3,h=6.25,k=1.5,plot=T)
summary(res3)           
res3$violations
x_13 <- PCA14[res3$violations,]

PCA15 <- PCA14[-res3$violations,]
dta_cusum.mqcd4 <- mqcd(PCA15)
res4 <- mqcs.mcusum(dta_cusum.mqcd4,h=6.25,k=1.5,plot=T)
summary(res4)           
res4$violations


nrow(PCA15)

## concatenate all removed points ##
removed <- rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x_10,x_11,x_12,x_13)
dim(removed)

match(x1[1:14,],PCA)
match(x2,PCA)
match(x3,PCA)
match(x4,PCA)
match(x5,PCA)
match(x6,PCA)
match(x7,PCA)
match(x8,PCA)
match(x9,PCA)
match(x_10,PCA)
match(x_11,PCA)
match(x_12,PCA)
match(x_13,PCA)
red_data <- data[-c(6,269,270,445,452,530,531,532,533,534,535,536,537,538,5,209,458,539,4,7,8,9,10,11,12,13,14,15,16,17,18,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,51,53,54,55,56,57,334,3,479,502,503,504,478,501,79,477,447,455,19,453,456,449,450,446,451,448,454),]
write.csv(red_data,"mean.csv")
Mean <- as.matrix(colMeans(red_data))
sigma_614 <- as.matrix(var(red_data))
write.csv(sigma_614,"sigma.csv")

### 462 points are retained as in control if we use t2 chart and then mcusum###
# NO points violated in mcusum, if applied after t2 iterations####

#######Mcusum##########(tried to code it manually)
k <- 1.5


mc_i[i] <- max(0,(sqrt(t(c_i)%*% solve(sigma123)%*% c_i)- (k*n_i[i+1])))


if (mc_i>0){
n_i[i+1] <- n_i[i]+1
} else {
  n_i[i+1] <- 1 }



### Mcusum using qcr #####
library(qcr)


dta_cusum.mqcd <- mqcd(PCA)
res <- mqcs.mcusum(dta_cusum.mqcd,h=7.5,k=2,plot=T)
summary(res)           
res$violations

PCA12 <- PCA[-res$violations,]
dta_cusum.mqcd1 <- mqcd(PCA12)
res1 <- mqcs.mcusum(dta_cusum.mqcd1,h=7.5,k=2,plot=T)
summary(res1)           
res1$violations


PCA13 <- PCA12[-res1$violations,]
dta_cusum.mqcd2 <- mqcd(PCA13)
res2 <- mqcs.mcusum(dta_cusum.mqcd2,h=7.5,k=2,plot=T)
summary(res2)           
res2$violations

PCA14 <- PCA13[-res2$violations,]
dta_cusum.mqcd3 <- mqcd(PCA14)
res3 <- mqcs.mcusum(dta_cusum.mqcd3,h=7.5,k=2,plot=T)
summary(res3)           
res3$violations

###### use below only when changing k #########

PCA15 <- PCA14[-res3$violations,]
dta_cusum.mqcd4 <- mqcd(PCA15)
res4 <- mqcs.mcusum(dta_cusum.mqcd4,h=7.5,k=2,plot=T)
summary(res4)           
res4$violations

PCA16 <- PCA15[-res4$violations,]
dta_cusum.mqcd5 <- mqcd(PCA16)
res5 <- mqcs.mcusum(dta_cusum.mqcd5,h=7.5,k=2,plot=T)
summary(res5)           
res5$violations




### T2 after cusum ######

PCA14
mu14 <- colMeans(PCA14)
var_pca14 <- var(PCA14)
dist14 <- c(rep(0,nrow(PCA14)))

for (i in 1:nrow(PCA14)){
  dist14[i] <- t(PCA14[i,]-mu14) %*% solve(var_pca14) %*% (PCA14[i,]-mu14)
}

dist14[dist14 > 11.07]
sum(dist14 > 11.07)

###2nd Iteration #####

PCAt2<- PCA14[dist14 < 11.07,]
mut2 <- colMeans(PCAt2)
var_pcat2 <- var(PCAt2)
distt2 <- c(rep(0,nrow(PCAt2)))

for (i in 1:nrow(PCAt2)){
  distt2[i] <- t(PCAt2[i,]-mut2) %*% solve(var_pcat2) %*% (PCAt2[i,]-mut2)
}

distt2[distt2 > 11.07]
sum(distt2 > 11.07)

### 3rd iteration ######

PCAt3<- PCAt2[distt2 < 11.07,]
mut3 <- colMeans(PCAt3)
var_pcat3 <- var(PCAt3)
distt3 <- c(rep(0,nrow(PCAt3)))

for (i in 1:nrow(PCAt3)){
  distt3[i] <- t(PCAt3[i,]-mut3) %*% solve(var_pcat3) %*% (PCAt3[i,]-mut3)
}

distt3[distt3 > 11.07]
sum(distt3 > 11.07)

## 4 th iteration

PCAt4<- PCAt3[distt3 < 11.07,]
mut4 <- colMeans(PCAt4)
var_pcat4<- var(PCAt4)
distt4 <- c(rep(0,nrow(PCAt4)))

for (i in 1:nrow(PCAt4)){
  distt4[i] <- t(PCAt4[i,]-mut4) %*% solve(var_pcat4) %*% (PCAt4[i,]-mut4)
}

distt4[distt4 > 11.07]
sum(distt4 > 11.07)

## 5th iteration ##

PCAt5<- PCAt4[distt4 < 11.07,]
mut5 <- colMeans(PCAt5)
var_pcat5<- var(PCAt5)
distt5 <- c(rep(0,nrow(PCAt5)))

for (i in 1:nrow(PCAt5)){
  distt5[i] <- t(PCAt5[i,]-mut5) %*% solve(var_pcat5) %*% (PCAt5[i,]-mut5)
}

distt5[distt5 > 11.07]
sum(distt5> 11.07)

### 6th iteration ####

PCAt6<- PCAt5[distt5 < 11.07,]
mut6 <- colMeans(PCAt6)
var_pcat6<- var(PCAt6)
distt6 <- c(rep(0,nrow(PCAt6)))

for (i in 1:nrow(PCAt6)){
  distt6[i] <- t(PCAt6[i,]-mut6) %*% solve(var_pcat6) %*% (PCAt6[i,]-mut6)
}

distt6[distt6 > 11.07]
sum(distt6> 11.07)

###7th iteration ####

PCAt7<- PCAt6[distt6 < 11.07,]
mut7 <- colMeans(PCAt7)
var_pcat7<- var(PCAt7)
distt7 <- c(rep(0,nrow(PCAt7)))

for (i in 1:nrow(PCAt7)){
  distt7[i] <- t(PCAt7[i,]-mut7) %*% solve(var_pcat7) %*% (PCAt7[i,]-mut7)
}

distt7[distt7 > 11.07]
sum(distt7> 11.07)

###8th iteration ####

PCAt8<- PCAt7[distt7 < 11.07,]
mut8 <- colMeans(PCAt8)
var_pcat8<- var(PCAt8)
distt8 <- c(rep(0,nrow(PCAt8)))

for (i in 1:nrow(PCAt8)){
  distt8[i] <- t(PCAt8[i,]-mut8) %*% solve(var_pcat8) %*% (PCAt8[i,]-mut8)
}

distt8[distt8 > 11.07]
sum(distt8> 11.07)

###### 395 variables are in control if we use t2 after mcusum ########


############## END #################

#### Phase 1 parameters from multi t2 #####
mu8
var_pca8

####Plot###
plot(dist8,type="b")
abline(h=11.07,col="red",lwd=3,lty=2)

library(IQCC)
cchart.T2.1(as.matrix(dist),552,1,5)

PCA1 <- PCA[dist < 11.07,]


### uni variate x bar charts after iterations ###

plot(PCA16[,1],type="b",main="Univariate X bar chart for PCA1 (Final Iteration)",ylab="PCA1",ylim=c(-6800,-6100))
abline(h= mean(PCA16[,1]),col="blue",lwd=3,lty=2)
abline(h= mean(PCA16[,1]) + (qnorm(1-.0009/2)*sqrt(var(PCA16[,1]))),col="red",lwd=3,lty=2)
abline(h= mean(PCA16[,1]) - (qnorm(1-.0009/2)*sqrt(var(PCA16[,1]))),col="red",lwd=3,lty=2)

plot(PCA16[,2],type="b",main="Univariate X bar chart for PCA2 (Final Iteration)",ylab="PCA2",ylim=c(30,500))
abline(h= mean(PCA16[,2]),col="blue",lwd=3,lty=2)
abline(h= mean(PCA16[,2]) + (qnorm(1-.0009/2)*sqrt(var(PCA16[,2]))),col="red",lwd=3,lty=2)
abline(h= mean(PCA16[,2]) - (qnorm(1-.0009/2)*sqrt(var(PCA16[,2]))),col="red",lwd=3,lty=2)

plot(PCA16[,3],type="b",main="Univariate X bar chart for PCA3 (Final Iteration)",ylab="PCA3",ylim=c(500,800))
abline(h= mean(PCA16[,3]),col="blue",lwd=3,lty=2)
abline(h= mean(PCA16[,3]) + (qnorm(1-.0009/2)*sqrt(var(PCA16[,3]))),col="red",lwd=3,lty=2)
abline(h= mean(PCA16[,3]) - (qnorm(1-.0009/2)*sqrt(var(PCA16[,3]))),col="red",lwd=3,lty=2)



