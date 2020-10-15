########################################################################################
############################## MORPHOLOGICAL VARIATION #################################

morpho_data<- read.table("morpho_total.csv", header = T, sep= "," )
data_raw<-morpho_data[,22:37] ## raw data
data_log<- morpho_data[,6:21] ##logarithmic data
species<- morpho_data[ ,5] 

############################ Discriminant Function Analysis ############################


#install.packages("MASS")
library(MASS)

LDA_table<- data.frame(species, data_log)

fit_LDA <- lda(LDA_table[,1] ~ .,data= LDA_table[,2:17])

loadings_LDA<-(fit_LDA$scaling)

scores_LDA <-predict(fit_LDA)$x


tiff("LDA.tiff", width=20, height=15, unit="cm", res=300)
quartz.options(height=10, width=12, dpi=72);
plot.new();
par(oma=c(1,1,1,1));
par(mar=c(3,3,1,1));
plot.window(xlim=c(-7,5), ylim=c(-4,4));

points(scores_LDA[110:262,1],scores_LDA[110:262,2],  col = '#FF6600', cex = 1.5, pch=16) #G
points(scores_LDA[263:436,1],scores_LDA[263:436,2],  col = '#D4AA00', cex = 1.5, pch=16) #H
points(scores_LDA[437:541,1],scores_LDA[437:541,2],  col = '#2CA02C', cex = 1.5, pch=16) #I
points(scores_LDA[101:109,1],scores_LDA[101:109,2],  col = '#7F2AFF', cex = 1.5, pch=16) #F
points(scores_LDA[85:100,1],scores_LDA[85:100,2],  col = '#800000', cex = 1.5, pch=16) #BCDE
points(scores_LDA[1:84,1],scores_LDA[1:84,2],  col = '#000080', cex = 1.5, pch=16) #A

axis(1, at=seq(-7, 5, by=3), cex.axis=1);
axis(2, at=seq(-4, 4, by=2), cex.axis=1, las=1);

mtext(side=1, text='DF1(60.05%)',line=2.5, cex=1)
mtext(side=2, text='DF2(24.21%)', line=2.8, cex=1)
legend<- c('A','BCDE','F','G','H','I') 
legend("topleft", legend = legend,y.intersp=0.8, ncol=2, cex=0.7, bty="n", col= c('#000080',
'#800000','#7F2AFF','#FF6600','#D4AA00','#2CA02C'), pch= 16)
dev.off()

write.csv(scores_LDA, "scores_LDA.csv") #saving LDA scores



########################################################################################


############################## Descriptive Statistics ##################################


Descriptive_table<- data.frame(species, data_raw)

A<- Descriptive_table[1:84,]
BCDE<- Descriptive_table[85:100,]
F <- Descriptive_table[101:109,]
G <- Descriptive_table[110:262,]
H <- Descriptive_table[263:436,]
I <- Descriptive_table[437:541,]


##  Media of each variable

media_A<- apply(A[,-1], 2, mean)
media_BCDE<- apply(BCDE[,-1], 2, mean)
media_F<- apply(F[,-1], 2, mean)
media_G<- apply(G[,-1], 2, mean)
media_H<- apply(H[,-1], 2, mean)
media_I<- apply(I[,-1], 2, mean)

##  Standard error of each variable

standard_error<- function(x) {
  2*(sd(x)/sqrt(length(x)))
}

A_error<-apply(A[,-1], 2, standard_error)
BCDE_error<-apply(BCDE[,-1], 2, standard_error)
F_error<-apply(F[,-1], 2, standard_error)
G_error<-apply(G[,-1], 2, standard_error)
H_error<-apply(H[,-1], 2, standard_error)
I_error<-apply(I[,-1], 2, standard_error)


A_min<- apply(A[,-1], 2, min)
BCDE_min<- apply(BCDE[,-1], 2, min)
F_min<- apply(F[,-1], 2, min)
G_min<- apply(G[,-1], 2, min)
H_min<- apply(H[,-1], 2, min)
I_min<- apply(I[,-1], 2, min)

A_max<- apply(A[,-1], 2, max)
BCDE_max<- apply(BCDE[,-1], 2, max)
F_max<- apply(F[,-1], 2, max)
G_max<- apply(G[,-1], 2, max)
H_max<- apply(H[,-1], 2, max)
I_max<- apply(I[,-1], 2, max)


### BOXPLOT with the most important variables from LDA analysis

#install.packages("ggplot2")

library(ggplot2)

table_boxplot<- data.frame(species, morpho_data[,c(27,28,33)])

fill <- c('#000080',  '#800000', '#7F2AFF', '#2CA02C', '#D4AA00','#FF6600')
          

BP<-ggplot(plan_box, aes(x=species, y=BP)) +
  geom_boxplot(position=position_dodge(1), fill= fill, alpha = 0.5)+
  theme_classic()
BP

ggsave("bp",  plot =BP, device = "pdf")

BOC<-ggplot(plan_box, aes(x=species, y=BOC)) +
  geom_boxplot(position=position_dodge(1), fill= fill, alpha = 0.5)+
  theme_classic()
BOC

ggsave("BOC",  plot =BOC, device = "pdf")

LN<-ggplot(plan_box, aes(x=species, y=LN)) +
  geom_boxplot(position=position_dodge(1), fill= fill, alpha = 0.5)+
  theme_classic()
LN

ggsave("LN",  plot =LN, device = "pdf")
