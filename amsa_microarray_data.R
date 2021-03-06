### Kristen Michelle Nader
### r0771801
### kristenmichelle.nader@student.kuleuven.be


### COMMENTS: Some variables need to be removed from the enviroment to spead up the analysis
### variables such as series61  are very large. In addition, I am only using the 1st 50 variables
### 
install.packages("janitor")
install.packages("class")
install.packages("corrplot")
install.packages("tree")

library("janitor")
library("class")
library("ggplot2")
library("gridExtra")
library("corrplot")
library("GEOquery")
library("grr")
library("leaps")
library("tree")
library("ROCR")
library("biomaRt")
library("cluster")


## logistic regression
data=read.table("C:\\Users\\user\\Desktop\\amsa_data\\GSE63061_reduced2.txt",sep=" ",header=TRUE)

### change gene names
series61 <- getGEO(filename = "C:\\Users\\user\\Desktop\\amsa_data\\GSE63061_series_matrix.txt")
GSE63061_expression <- as.data.frame(t(assayData(series61)$exprs))

probes <- colnames(GSE63061_expression)


probeData61 <- as.data.frame(matrix(nrow = length(probes)))
probeData61$probes <- probes
head(probeData61$probes)
probeData61$EntrezID <- featureData(series61)[["Entrez_Gene_ID"]][matches(probeData61$probes,featureData(series61)[[1]], all.y = F)$y]
probeData61$GeneName <- featureData(series61)[["Symbol"]][matches(probeData61$probes,featureData(series61)[[1]], all.y = F)$y]


k=colnames(data)[-1]
l=probeData61$probes
for (val in 1:length(k)) {
  for (j in 1:length(probeData61$probes))
    if (k[val]==probeData61$probes[j])
      # print(k[val])
      k[val]=probeData61$GeneName[j]
}

colnames(data)=c("class",k)
data <- clean_names(data)

## clean up workspace series 61 is a large expression set and will freeze the computer( froze my comp) 
rm(series61,GSE63061_expression,probes,probeData61,l,val,j)

data=data[,c(1:50)]


## transform to binary and use logistic regression


data$class=ifelse(data$class=="CTL", "CTL","EXP")



#training=Data[c(1:1000),c(1:50)]
#testing=Data[c(1001:1500),c(1:50)]
training=data[c(1:254),c(1:50)]
testing=data[c(255:382),c(1:50)]
cutoff=0.5
set.seed(1)





glm1=glm(training$class~.,family=binomial,data=training)
summary(glm1)

slm.back<- step(glm(as.factor(training$class) ~., family=binomial,training), direction = "backward")
summary(slm.back)


glm1_data_pred <- predict(slm.back, testing, type = "response")
predicted.glm1 <- as.numeric(glm1_data_pred > cutoff)
table(testing$class, predicted.glm1)
missclassification_logistic=(22+23)/(15+22+23+68)

predict<-fitted(slm.back)
pred<-prediction(predict,training$class)
perf<-performance(pred,measure="tpr",x.measure="fpr")
plot(perf,main="sensitivity vs false positive rate",colorize=TRUE,
     colorkey.relwidth=0.5,lwd=4.5)
perf_auc=performance(pred,measure="auc")
perf_auc@y.values



## using the set of identifiers, create a tree looking for structure in ctl and exp


data.tree1 <- tree(as.factor(training$class)~hla_a29_1+rpl26+loc388588+folr3+ifi44l
                   +loc648622+cox7c+loc441377
                   +loc653658+ccdc72+hbe1+rpl23+loc644934
                   +rps26+hint1+klrb1
                     ,method="recursive.partition",split="deviance",data=training)

tree.control(nobs=254)
attributes(data.tree1)
summary(data.tree1)
data.tree1
plot(data.tree1) 
text(data.tree1,splits=T,all=T,pretty = 0)              # puts names on tree plot
title("Tree selection of variables-Backward selection variables ")

## check misclass of this tree
tree1.pred=predict(data.tree1,testing,type="class")
table(tree1.pred,testing$class)
missclass_tree1=(20+29)/(17+20+29+62)

## misclassification of the pruned tree

data.tree.cv1 <- cv.tree(data.tree1,FUN=prune.tree)
plot(data.tree.cv1$size,data.tree.cv1$k,type="l")
plot(data.tree.cv1$size,data.tree.cv1$dev,type="l",xlab="Number of end nodes",ylab="Deviance") 
T.ind1<-prune.tree(data.tree1,best=5) 
summary(T.ind1)

plot(T.ind1)
text(T.ind1, col="purple")
title("Pruned Tree using CTL/EXP and Backward selection variables- 5 terminal nodes") 

## missclass
T.ind1.pred=predict(T.ind1,testing,type="class")
table(T.ind1.pred,testing$class)
missclass_tree_T.ind1=(31)/(7+1+30+90)



T.ind10<-prune.tree(data.tree1,best=10) 
summary(T.ind10)

plot(T.ind10)
text(T.ind10,col="purple")
title("Pruned Tree using CTL/EXP and Backward selection variables- 10 terminal nodes") 


T.ind10.pred=predict(T.ind10,testing,type="class")
table(T.ind10.pred,testing$class)
missclass_tree_T.ind10=(7+26)/(11+7+26+84)


## try using the full data variables with the training/testing
data.tree2 <- tree(as.factor(training$class) ~ ., data=training, method="recusive.partition",
                   split="deviance")
tree.control(nobs=254)

attributes(data.tree2)
summary(data.tree2)

data.tree2
plot(data.tree2)
text(data.tree2,splits=T,all=T)              # puts names on tree plot
title("tree selection of variables-full number of variables(50)")

## misclassif

tree2.pred=predict(data.tree2,testing,type="class")
table(tree2.pred,testing$class)
missclass_tree_2=(28+17)/(20+17+28+63)




data.tree.cv2 <- cv.tree(data.tree2,FUN=prune.tree)
plot(data.tree.cv2$size,data.tree.cv2$k,type="l")
plot(data.tree.cv2$size,data.tree.cv2$dev,type="l",xlab="Number of end nodes",ylab="Deviance")

T.ind2<-prune.tree(data.tree2,best=5)

summary(T.ind2)

par(mfrow=c(1,1))
plot(T.ind2)
text(T.ind2,pretty = 1,col="purple")
title("Pruned tree-full variables(50)-5 terminal") ## now the error is 0.21 but here we tried using the entire dataset

T.ind2.pred=predict(T.ind2,testing,type="class")
table(T.ind2.pred,testing$class)
missclass_tree_2=(30+0)/(7+0+30+91)

T.ind3<-prune.tree(data.tree2,best=10)

summary(T.ind3)

par(mfrow=c(1,1))
plot(T.ind3)
text(T.ind3,pretty = 1,col="purple")
title("Pruned tree-full variables(50)-10 terminal") ## now the error is 0.21 but here we tried using the entire dataset

T.ind3.pred=predict(T.ind3,testing,type="class")
table(T.ind3.pred,testing$class)
missclass_tree_3=(4+26)/(11+4+26+87)







### use clustering 

# hla_a29_1+rpl26+loc388588+folr3+ifi44l
# +loc648622+cox7c+loc441377
# +loc653658+ccdc72+hbe1+rpl23+loc644934
# +rps26+hint1+klrb1

data_subset=training[,c(3,10,13,16,18,20,22,25,27,30,31,37,40,44,45,50)]

data.clusD=hclust(dist(data_subset),method="ward.D")
plot(data.clusD,xlab="Data",ylab="Ward.D methode",sub="")



plot(c(1,15),c(-0,20),type="n",xlab="var",ylab="Value",main="Profile Plot - Expression Data")
# Use a loop to generate profiles for each observation.
for (k in (1:254))
{
  if(training$class[k]=="EXP")
    points(1:15,data_subset[k,],type="l",col=1)
  else
    points(1:15,data_subset[k,],type="l",col=2)
}

data_subset_cluster <- cutree(data.clusD,k=2)
table(data_subset_cluster)
table(training$class)


pairs(data_subset,col=data_subset_cluster,main="Predicted grouping structure in the data")  
pairs(data_subset,col=as.factor(data$Class),main="Observed grouping structure in the data")  



p1 <- ggplot(data, aes(x=hla_a29_1, y=rpl26, col=as.factor(data$class))) +
  geom_point(alpha=0.8)

p2 <- ggplot(data, aes(x=loc653658, y=snrpg, col=as.factor(data$class))) +
  geom_point(alpha=0.8)

p3 <- ggplot(data, aes(x=loc653658, y=ccdc72, col=as.factor(data$class))) +
  geom_point(alpha=0.8)

grid.arrange(p1 ,ncol=1)
print(p1 + ggtitle("Expression of hla_a29_1(x) and rpl26(y)"))

grid.arrange(p2 ,ncol=1)
print(p2 + ggtitle("Expression of loc653658(x) and snrpg(y)"))

grid.arrange(p3 ,ncol=1)
print(p3 + ggtitle("Expression of loc653658(x) and ccdc72(y)"))
          

## KNN

train_temp=training[,-1]
test_temp=testing[,-1]

train_acc <- c()


for (i in 1:5){
  set.seed(1)
  train_pred <- knn(train_temp, train_temp, training$class, k=i)
  train_acc <- c(train_acc, mean(train_pred == training$class))
}

tep_df<- data.frame(NN=c(1:5), train_accu=train_acc)

  
p4<- ggplot(data=tep_df, aes(x=NN, y=train_accu, group=1)) +
  geom_line(color="pink")+
  geom_point()


grid.arrange(p4 ,ncol=1)
print(p4 + ggtitle("Training Accuracy for k-NN"))


train_pred <- knn(train_temp, train_temp, training$class, k=3)
accuracy <- mean(train_pred == training$class)
cat("Training Accuracy: ", accuracy, sep='')

testing_pred=knn(train_temp, test_temp, training$class, k=3)
testing_accuracy=mean(testing_pred==testing$class)
cat("Testing Accuracy: ", testing_accuracy, sep='')


ensembl_human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
attrib<-listAttributes(ensembl_human)



structureData_human <- getBM(mart = ensembl_human,
                             attributes = c("definition_1006"),
                             filters = "external_gene_name",
                             values = c("LOC653658"))
head(structureData_human)



