# Packages used GGally,ggplot2,gridExtra,car,carData,MVN,heplots,matrixcalc,psych,factoextra,FactoMineR,MASS
# Data is stored in the file "~//Stock Data Assignment//stocks.txt" relative to the present
# working directory (can be obtained using getwd() command in R console)
# So if getwd() gives output "C:/Users/Soumya/Documents"
# Then full file path of the data is "I:/bankruptcy.txt"
# Please change the file location as required

# Please remove the Hash sign from the code lines with install.packages
# if the packages are not installed 

# Installing/loading required packages
#install.packages(c("GGally","ggplot2","gridExtra","car","carData","MVN","heplots","matrixcalc","psych","factoextra","FactoMineR","MASS"))
library(GGally)
library(ggplot2)
library(gridExtra)
library(car)
library(carData)
library(MVN)
library(heplots)
library(matrixcalc)
library(psych)
library(factoextra)
library(FactoMineR)
library(MASS)


## Importing the data
#setwd("C:/Users/SUSMIT/Desktop")
#bankruptcy=read.csv(file="bankruptcy.csv")
bankruptcy=read.delim("~/bankruptcy.txt",sep="",header=FALSE)


## Assigning variables to columns of the data frame
data=data.frame(bankruptcy)
b = split(data, data[,5])
br = as.matrix(b$`0`[,1:4]) # Bankrupt firms
fs = as.matrix(b$`1`[,1:4]) # Financially sound firms
n = c(nrow(br), nrow(fs))
financialcond=as.factor(data[,5])


## Exploratory Data Analysis(original dataset)
ggpairs(data, columns =1:4,columnLabels = c("X1","X2","X3","X4"), ggplot2::aes(colour=as.factor(data[,5])))
p1=ggplot(data, aes(y=data[,1],fill=financialcond))+
  geom_boxplot()+labs(title="boxplot for X1",y="X1",fill="Financial\nCondition")
p2=ggplot(data, aes(y=data[,2],fill=financialcond))+
  geom_boxplot()+labs(title="boxplot for X2",y="X2",fill="Financial\nCondition")
p3=ggplot(data, aes(y=data[,3],fill=financialcond))+
  geom_boxplot()+labs(title="boxplot for X3",y="X3",fill="Financial\nCondition")
p4=ggplot(data, aes(y=data[,4],fill=financialcond))+
  geom_boxplot()+labs(title="boxplot for X4",y="X4",fill="Financial\nCondition")
grid.arrange(p1,p2,p3,p4,nrow=2,ncol=2)



## PCA

## Bankrupt group
pca_br=PCA(br, scale.unit = FALSE, ncp = 5, graph = FALSE)
eig_br=pca_br$eig
fviz_eig(pca_br, addlabels = TRUE, ylim = c(0, 75))
prcomp(br,scale. = FALSE)

## Financially sound group
pca_fs=PCA(fs, scale.unit = FALSE, ncp = 5, graph = FALSE)
eig_fs=pca_fs$eig
fviz_eig(pca_fs, addlabels = TRUE, ylim = c(0, 95))
prcomp(fs,scale. = FALSE)

## Checking for multivariate normality
mvn(br,mvnTest = "royston",multivariatePlot = "qq")
mvn(fs,mvnTest = "royston",multivariatePlot = "qq")


## Power transformation
powerTransform(fs,family = "bcnPower")
fs_trans=bcnPower(fs,lambda = c( 0.43351945,1.01352799,0.04555263,1.99833949),gamma = c(0.2359694,0.1664129,3.7527159,3.0851024))
mvn(fs_trans,mvnTest = "royston",multivariatePlot = "qq")
br_trans=bcnPower(br,lambda = c( 0.43351945,1.01352799,0.04555263,1.99833949),gamma = c(0.2359694,0.1664129,3.7527159,3.0851024))


## Checking multivariate normality of the transformed data
mvn(br_trans,mvnTest = "royston")
mvn(fs_trans,mvnTest = "royston")


## Transformed data matrix
d1=rbind(br_trans,fs_trans)
c=rep(c(0,1),n)
data_trans=cbind(d1,c)

## Exploratory data analysis of the transformed data
ggpairs(data.frame(data_trans), columns = 1:4,columnLabels = c("X1","X2","X3","X4"), ggplot2::aes(colour=as.factor(data[,5])))
p21=ggplot(data.frame(data_trans), aes(y=data_trans[,1],fill=financialcond))+
  geom_boxplot()+labs(title="boxplot for X1",y="X1",fill="Financial\nCondition")
p22=ggplot(data.frame(data_trans), aes(y=data_trans[,2],fill=financialcond))+
  geom_boxplot()+labs(title="boxplot for X2",y="X2",fill="Financial\nCondition")
p23=ggplot(data.frame(data_trans), aes(y=data_trans[,3],fill=financialcond))+
  geom_boxplot()+labs(title="boxplot for X3",y="X3",fill="Financial\nCondition")
p24=ggplot(data.frame(data_trans), aes(y=data_trans[,4],fill=financialcond))+
  geom_boxplot()+labs(title="boxplot for X4",y="X4",fill="Financial\nCondition")
grid.arrange(p21,p22,p23,p24,nrow=2,ncol=2)

## Check the equality of the dispersion matrices
boxM(data_trans[,-5],as.factor(data_trans[,5]))


## Check for the equality of the mean vectors
m1_trans = apply(br_trans, 2, mean)
m2_trans = apply(fs_trans, 2, mean) 
s1_trans=var(br_trans)
s2_trans=var(fs_trans)
T_sq=t(m1_trans-m2_trans)%*%solve(s1_trans/n[1]+s2_trans/n[2])%*%(m1_trans-m2_trans)
m1=matrix.power(((s1_trans/n[1])%*%solve(s1_trans/n[1]+s2_trans/n[2])),2)
m2=matrix.power(((s2_trans/n[2])%*%solve(s1_trans/n[1]+s2_trans/n[2])),2)
m3=(s1_trans/n[1])%*%solve((s1_trans/n[1])+(s2_trans/n[2]))
m4=((s2_trans/n[2])%*%solve(s1_trans/n[1]+s2_trans/n[2]))
v_nu=(tr(m1)+tr(m3)^2)/n[1]+(tr(m2)+tr(m4)^2)/n[2]
v=20/v_nu
crit=(43*4/(43-4+1))*qf(0.95,4,43-4+1)


## Interval estimation of mean differences

## Separate interval
dml = (m1_trans-m2_trans) - qt(0.975,n[1]+n[2]-2)*sqrt(((1/n[1])+(1/n[2]))*((n[1]-1)*diag(s1_trans) + (n[2]-1)*diag(s2_trans))/(n[1]+n[2]-2))
dmu=(m1_trans-m2_trans)+ qt(0.975,n[1]+n[2]-2)*sqrt(((1/n[1])+(1/n[2]))*((n[1]-1)*diag(s1_trans) + (n[2]-1)*diag(s2_trans))/(n[1]+n[2]-2))                                                       
cbind(dml,dmu)
abs(dml-dmu)

## Simultaneous interval
dml_1=(m1_trans-m2_trans)-sqrt(qchisq(0.95,4))*sqrt(diag(s1_trans)/n[1]+diag(s2_trans)/n[2])
dmu_1=(m1_trans-m2_trans)+sqrt(qchisq(0.95,4))*sqrt(diag(s1_trans)/n[1]+diag(s2_trans)/n[2])
cbind(dml_1,dmu_1)
abs(dml_1-dmu_1)

## Factor analysis:

## Bankrupt group

pca_br_trans=PCA(br_trans, scale.unit = FALSE, ncp = 5, graph = FALSE)
eig_br_trans=pca_br_trans$eig
fviz_eig(pca_br_trans, addlabels = TRUE, ylim = c(0, 75))
prcomp(br_trans,scale. = FALSE)

fa_br_trans=fa(br_trans,nfactors = 2,rotate = 'none',fm='pa')
fa_br_trans

fa_bartlett_br_trans=fa(br_trans,nfactors = 2,rotate = 'none',fm='pa',scores="Bartlett")
scores_bartlett_br_trans=fa_bartlett_br_trans$scores
mvn(scores_bartlett_br_trans,mvnTest="royston")
cor.test(scores_bartlett_br_trans[,1],scores_bartlett_br_trans[,2])

fa_regression_br_trans=fa(br_trans,nfactors = 2,rotate = 'none',fm='pa',scores="regression")
scores_regression_br_trans=fa_regression_br_trans$scores
mvn(scores_regression_br_trans,mvnTest="royston")
cor.test(scores_regression_br_trans[,1],scores_regression_br_trans[,2])



## Financially sound group

pca_fs_trans=PCA(fs_trans, scale.unit = FALSE, ncp = 5, graph = FALSE)
eig_fs_trans=pca_fs_trans$eig
fviz_eig(pca_fs_trans, addlabels = TRUE, ylim = c(0, 95))
prcomp(fs_trans,scale. = FALSE)


fa_fs_trans=fa(fs_trans,nfactors = 2,rotate = 'none',fm='pa')
fa_fs_trans


fa_bartlett_fs_trans=fa(fs_trans,nfactors = 2,rotate = 'none',fm='pa',scores="Bartlett")
scores_bartlett_fs_trans=fa_bartlett_fs_trans$scores
mvn(scores_bartlett_fs_trans,mvnTest="royston")
cor.test(scores_bartlett_fs_trans[,1],scores_bartlett_fs_trans[,2])

fa_regression_fs_trans=fa(fs_trans,nfactors = 2,rotate = 'none',fm='pa',scores="regression")
scores_regression_fs_trans=fa_regression_fs_trans$scores
mvn(scores_regression_fs_trans,mvnTest="royston")
cor.test(scores_regression_fs_trans[,1],scores_regression_fs_trans[,2])



## Discriminant analysis
qda <- qda(data_trans[,-5],data_trans[,5],prior=c(0.5,0.5))
qda.pred <- predict(qda)$class
table(data_trans[,5], qda.pred, dnn = c('Actual Group','Predicted Group'))
aper=(4*100)/46


qda.cv <- qda(data_trans[,-5],data_trans[,5],prior=c(0.5,0.5), CV = TRUE)
table(as.factor(data_trans[,5]), qda.cv$class, dnn = c('Actual Group','Predicted Group'))
holdoutaper=(8*100)/46




