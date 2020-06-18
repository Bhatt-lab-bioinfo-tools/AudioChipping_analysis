#remove NEB from fig3 and 4, generate fig 3 1 (AT vs CT, all phenes), and fig 3 2 (non NA phenes AT,CT versus NT)
set.seed(123)
do_imputation=FALSE

#import spss data
library(foreign)
#boxcox
library(car)
library(MASS)

library(bestNormalize)


library(raster)

myClip <- function(x, a=-2, b=2) {
a<-ifelse(x <= a,  a, ifelse(x >= b, b, x))
return(as.numeric(a))
}


library(ggdendro)
# Do something, or tell me why it failed



my_update_function <- function(x){
    tryCatch(
        # This is what I want to do...
        {
        y = x * 2
        return(y)
        },
        # ... but if an error occurs, tell me what happened: 
        error=function(error_message) {
            message("This is my custom message.")
            message("And below is the error message from R:")
            message(error_message)
            return(NA)
        }
    )
}

#if lsr library is not available, we can manually calculate cohen's D
cohens_d <- function(x, y) {
    lx <- length(x)- 1
    ly <- length(y)- 1
    md  <- (mean(x) - mean(y))        ## mean difference (numerator)
    #md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
    csd <- lx * var(x) + ly * var(y)
    csd <- csd/(lx + ly)
    csd <- sqrt(csd)                     ## common sd computation

    cd  <- md/csd                        ## cohen's d
}
#define Min-Max normalization function
normalize <- function(x) {
    return ((x - min(x)) / (max(x) - min(x)))
}

#apply Min-Max normalization to first four columns in iris dataset
normalize2 <- function(x) {
   x <- (x - min(x)) / (max(x) - min(x))
   #return (scale(x))
   return (scale(x, center=TRUE, scale=TRUE))
   #return ((x-mean(x))/sd(x))
}

#if lsr library is available, we can just load this library and directly calculate cohen's D
library(lsr) #cohen's D, effect size, cohensD(x,y)

#library for heatmap plot
library(gplots)
#bluered collor pallete
library(squash)
#create pallete from red-blue scale
colpal <- bluered(15)

#two-way clustering library
library(mixOmics)

#load dataset from spss sav file format
dataset = read.spss("mydata.sav", to.data.frame=TRUE)
print("BEFORE_impute")
nrow(dataset)

write.table(dataset, file="original_dataset.csv", sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)

#dataset <- na.omit(dataset)

library(mice)

if(do_imputation==FALSE){
    dataset <- na.omit(dataset)
    suffix="imputationFALSE"
}else{
    mice_imputes <- mice(dataset, m=5, maxit = 40)
    dataset <- complete(mice_imputes,5)
    suffix = "imputationTRUE"
}
raw <- dataset$SSQ12_average

dataset$SSQ12_average <- NULL

tmp <- dataset
print("AFTER_impute")
nrow(dataset)

#shapiro test for normality
dataset2 <- dataset
start=3
lshap <- lapply(as.data.frame(tmp[,start:ncol(tmp)]), shapiro.test)
lres <- sapply(lshap, `[`, c("statistic","p.value"))
t(lres)

lshap <- lapply(as.data.frame(dataset[,start:ncol(dataset)]), shapiro.test)
lres1 <- sapply(lshap, `[`, c("statistic","p.value"))

print("AFTER_COPY")
cbind(as.data.frame(t(lres)), as.data.frame(t(lres1)) )


#null hypothesis = data is normal
#p<0.05 reject null hypothesis
to_norm <- subset(as.data.frame(t(lres)), p.value<0.05)

to_norm <- row.names(to_norm)




#print(to_norm)
print("TO_NORMALIZE")
print(to_norm)

start=3
print(ncol(dataset))
for(i in c(start:ncol(dataset)) ){
   
    a <- 0
#<- scale(dataset[,i],center=TRUE,scale=TRUE)
    dataset[,i] <- scale(dataset[,i],center=TRUE,scale=TRUE)
    #dataset[,i] <- boxCoxVariable(dataset[,i])
    #print(min(dataset[,i]))

}
lshap <- lapply(as.data.frame(dataset[,start:ncol(dataset)]), shapiro.test)
lres1 <- sapply(lshap, `[`, c("statistic","p.value"))

print("AFTER_SCALING")
cbind(as.data.frame(t(lres)), as.data.frame(t(lres1)) )

fname=paste0("data_transform_", suffix, ".pdf")
pdf(file=fname)

par(mfrow = c(1,2))

for(i in to_norm){

    x <- dataset[,i]
    
    dataset[,i] <- as.numeric(dataset[,i]+abs(min(dataset[,i]))+0.001)
    #dataset[,i] <- boxcox(dataset[,i]~1, objective.name="Shapiro-Wilk")
    BNobject <- bestNormalize(dataset[,i])
    my_new_data <-  predict(BNobject, newdata = dataset[,i])
    dataset[,i] <- my_new_data
    
    xx <- seq(min(x), max(x), length = 100)
    MASS::truehist(x, 
                  main = paste0(i,": Before Transformation"), nbins = 10)
    MASS::truehist(BNobject$x.t, 
                  main = paste0("Best Transformation: ", 
                                class(BNobject$chosen_transform)[1]), nbins = 10)
    
    #plot(xx, predict(BNobject, newdata = xx), type = "l", col = 1, 
    #     main = "Best Normalizing transformation", ylab = "g(x)", xlab = "x")

    #dataset[,i] <- boxcox(dataset[,i]~1, objective.name="Log-Likelihood")
}

dev.off()


for(i in c(start:ncol(dataset)) ){
   
    a <- 0
#<- scale(dataset[,i],center=TRUE,scale=TRUE)
#    dataset[,i] <- scale(dataset[,i],center=TRUE,scale=TRUE)
    #dataset[,i] <- boxCoxVariable(dataset[,i])
    #print(min(dataset[,i]))

}



dataset2 <- dataset


print("AFTER NORMALIZATION AND BOX-COX TRANSFORMATION")
lshap <- lapply(as.data.frame(dataset[,start:ncol(dataset)]), shapiro.test)
lres2 <- sapply(lshap, `[`, c("statistic","p.value"))
cbind(as.data.frame(t(lres)),as.data.frame(t(lres2)))

#show a portion of the file, just for sanity check
print(dataset$Category)
summary(dataset)

fname=paste0("normalized_data_", suffix, ".csv")

write.table(dataset, file=fname, sep=",", quote=FALSE, row.names=FALSE, col.names=TRUE)

#split the dataset among the 3 groups
controls <- subset(dataset, Category=="Mid SSQ")
low_ssq <-  subset(dataset, Category=="Low SSQ")
high_ssq <-  subset(dataset, Category=="High SSQ")

mid_low <- subset(dataset, Category=="Low SSQ" | Category=="Mid SSQ")
mid_high <- subset(dataset, Category=="High SSQ" | Category=="Mid SSQ")

#These variables will store the effect size of cohen's D
control_vs_low_ssq_cohen <- NULL
control_vs_high_ssq_cohen <- NULL
high_ssq_vs_low_ssq_cohen <- NULL

j=1
#loop through all non ID;label columns and calculate cohen's D (effect size)

library("effsize")

a <- cohen.d(low_ssq[,start], controls[,start])
cohen.d(high_ssq[,start], controls[,start])
print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
print(a$estimate)

for(i in c(start:ncol(controls))){
    #print(cbind(controls[,i], low_ssq[,i]))    
    #control_vs_low_ssq_cohen[j] <- cohensD(controls[,i], low_ssq[,i], method="raw")
    #control_vs_high_ssq_cohen[j] <- cohensD(controls[,i], high_ssq[,i], method="raw")

    #control_vs_low_ssq_cohen[j] <- cohens_d(controls[,i], low_ssq[,i])
    #control_vs_high_ssq_cohen[j] <- cohens_d(controls[,i], high_ssq[,i])
    
    control_vs_low_ssq_cohen[j] <- cohen.d(low_ssq[,i], controls[,i])$estimate
    control_vs_high_ssq_cohen[j] <- cohen.d(high_ssq[,i], controls[,i])$estimate

    j=j+1
}


j=1

for(i in c(start:ncol(controls))){
    
    #high_ssq_vs_low_ssq_cohen[j] <- cohensD(high_ssq[,i], low_ssq[,i])
    high_ssq_vs_low_ssq_cohen[j] <- cohens_d(high_ssq[,i], low_ssq[,i])
   
    j=j+1
}



#Row clustering using Pearson correlation among phenes
hr <- hclust(as.dist(1-cor((dataset[,start:ncol(controls)]), method="pearson")), method="complete")
print("hr done")


#To cluster by column, we need more than 2 diseases!!! Leaving this commented by now
#hc <- hclust(as.dist(1-cor(XXXXXXXX, method="pearson")), method="complete")

#plot heatmap with dendogram
#data_to_plot <- as.matrix(data.frame(AT=control_vs_low_ssq_cohen, CT=control_vs_high_ssq_cohen))


data_to_plot <- as.matrix(data.frame(Low=control_vs_low_ssq_cohen, High=control_vs_high_ssq_cohen))

row.names(data_to_plot) <- colnames(dataset[,start:ncol(controls)])

fname=paste0("data_plotted_fig3_cohenD_", suffix, ".csv")

write.table(data_to_plot, file=fname, sep=",", quote=FALSE, row.names=TRUE, col.names=TRUE)

fname=paste0("fig3_normalized_", suffix, ".pdf")

pdf(file=fname, height=10, width=6)

heatmap.2(data_to_plot, Rowv=as.dendrogram(hr), col=bluered, scale="col", density.info="none", 
         trace="none", lhei = c(2.5, 10), margins=c(4,15), srtCol=0,  adjCol = c(0.5,1))

dev.off()

fname=paste0("fig3_normalized_rescaled_", suffix, ".pdf")

data_to_plot <- scale(data_to_plot)
for(i in 1:nrow(data_to_plot)){

    data_to_plot[i,] <- myClip( data_to_plot[i,], a=-2, b=2 )

}


pdf(file=fname, height=10, width=6)


heatmap.2(scale(data_to_plot), Rowv=as.dendrogram(hr), col=bluered, scale="none", density.info="none", 
         trace="none", lhei = c(2.5, 10), margins=c(4,15), srtCol=0,  adjCol = c(0.5,1))

dev.off()


fname=paste0("fig3_normalized_", suffix, ".png")

png(file=fname)

heatmap.2(scale(data_to_plot), Rowv=as.dendrogram(hr), col=bluered, scale="col", density.info="none", 
         trace="none", lhei = c(2.5, 10), margins=c(4,15), srtCol=0,  adjCol = c(0.5,1))

dev.off()


#Row clustering using Pearson correlation among phenes
hr <- hclust(as.dist(1-cor((dataset[,start:ncol(controls)]), method="pearson")), method="complete")
print("hr2 done")

#plot heatmap with dendogram
#data_to_plot <- as.matrix(data.frame(CT=high_ssq_vs_low_ssq_cohen, AT=high_ssq_vs_low_ssq_cohen))
data_to_plot <- as.matrix(data.frame(RT=high_ssq_vs_low_ssq_cohen, ST=high_ssq_vs_low_ssq_cohen))
print("data_to_plot figure 3")
data_to_plot
row.names(data_to_plot) <- colnames(dataset[,start:ncol(controls)])

#pdf(file="fig3_noNEB_CT-vs-AT_RD.pdf")

#fname=paste0("fig3_normalized_zscores_Pearson_correl_clustering_", suffix, ".pdf")

#pdf(file=fname)

#heatmap.2(data_to_plot, Rowv=as.dendrogram(hr), Colv=FALSE, col=bluered, scale="col", density.info="none", 
#         trace="none", lwid = c(2,4), lhei = c(2.5, 10), margins=c(8,20), srtCol=0,  adjCol = c(0.5,1), keysize=c(1))

#dev.off()

#fname=paste0("fig3_normalized_zscores_Pearson_correl_clustering_", suffix, ".png")

#heatmap.2(data_to_plot, Rowv=as.dendrogram(hr), Colv=FALSE, col=bluered, scale="col", density.info="none", 
#         trace="none", lwid = c(5, 8), lhei = c(2.5, 10), margins=c(8,20), srtCol=0,  adjCol = c(0.5,1), keysize=c(1))

#dev.off()


#convert data to matrix
start=3
#for(i in c(start:ncol(dataset2)) ){

 #  dataset2[,i] <- normalize2(dataset2[,i])


#}

#dataset2[,start:ncol(dataset2)] <- scale(dataset2[,start:ncol(dataset2)], scale=TRUE,center=TRUE)
#dataset2[,start:ncol(dataset2)] <- scale(dataset2[,start:ncol(dataset2)], scale=TRUE,center=TRUE)
#dataset2[,start:ncol(dataset2)] <- scale(dataset2[,start:ncol(dataset2)], scale=TRUE,center=TRUE)
dataset2 <- dataset
print(dataset2)

data_to_plot <- as.matrix(dataset2[,start:ncol(controls)])



#set row names (subjects)
row.names(data_to_plot) <- dataset[,1]

#set column names (phenes)
names(data_to_plot) <- names(dataset)[start:ncol(controls)]

#transpose
#data_to_plot <- scale(data_to_plot,center = TRUE, scale = TRUE)
data_to_plot = t(data_to_plot)
#data_to_plot <- scale(data_to_plot,center = TRUE, scale = TRUE)



#color subjects by NT, AT, CT
cond=as.character(dataset[,2])
cond_continuous<-as.numeric(as.character(raw))
cond_continuous_raw<-as.numeric(as.character(raw))

print(cond_continuous_raw[1:10])

cond_continuous<-scale(cond_continuous, center=TRUE, scale=FALSE)
cond_continuous <- myClip(cond_continuous, a=-1, b=1) 

print("CONTINUOUS")
print(cond_continuous)

print("SAMPLE COL BEFORE CLUSTERING")

print(cond_continuous_raw[1:10])

length(dataset[,1])
length(cond_continuous_raw)
length(cond_continuous)
length(cond)

print(data.frame(ID=dataset[,1],raw=cond_continuous_raw,condition=cond_continuous, category=cond))


cond.col <- c("Low SSQ" = "red", "Mid SSQ" = "black", "High SSQ" = "blue")
print("condition")
print(cond)
print("last_plot")

head(data_to_plot)

#create heatmap object
obj.cim=cim(data_to_plot)

#dissimilarity matrix
ddr <- as.hclust(obj.cim$ddr)

#clustering, cut tree in height level = 3
cl <- cutree(ddr, k = 3)

# Put the labels at the same height: hang = -1
plot(ddr, hang = -1, cex = 0.6)

fname=paste0("dendogram_", suffix, ".pdf")

p <- ggdendrogram(ddr, rotate = TRUE, size = 4, theme_dendro = FALSE)
pdf(fname, height=6,width=12)
plot(p)
dev.off()

a <- cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("correlation", "correlation"),
    margins = c(5, 10), keysize=c(1.5,1), color=colpal
   #lhei=c(1,0.1,1), scale=TRUE
    )


summary(a$mat)
str(a$mat)

fname=paste0("data_plotted_fig4_", suffix, ".csv")

write.table(a$mat, file=fname, quote=FALSE, sep=",")


#a$mat <- cor(data_to_plot)

nrow(data_to_plot)
nrow(t(data_to_plot))



lshap <- lapply(as.data.frame(t(data_to_plot)), shapiro.test)

lres <- sapply(lshap, `[`, c("statistic","p.value"))
t(lres)
print("original")

fname=paste0("fig4_normalized_zscores_euclidean_clustering_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("correlation", "correlation"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE
   #lhei=c(1,0.1,1), scale=TRUE
    )

dev.off()

fname=paste0("fig4_normalized_zscores_pearson_correl_clustering_", suffix, ".pdf")

pdf(file=fname, height=10,width=14)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("correlation", "correlation"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()



for(i in 1:nrow(data_to_plot)){

    data_to_plot[i,] <- myClip( data_to_plot[i,] )

}


fname=paste0("fig4_normalized_Zscores_cliped_euclidean_clustering_", suffix, ".pdf")


pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("correlation", "correlation"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()


fname=paste0("fig4_normalized_Zscores_cliped_correl_clustering_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("correlation", "correlation"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()

#the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of 
#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

fname=paste0("fig4_normalized_Zscores_cliped_correl_clustering_ward_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("correlation", "correlation"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE, clust.method = c("ward","ward")
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()


#kmeans(data_to_plot, centers, iter.max = 50, nstart = 1,
#       algorithm = c("Hartigan-Wong", "Lloyd", "Forgy",
#                     "MacQueen"), trace=FALSE)

#kc1_3 <- kmeans(data_to_plot, 3, iter.max = 50, nstart = 1,algorithm = c("Hartigan-Wong"), trace=FALSE)
#kc2_3 <- kmeans(data_to_plot, 3, iter.max = 50, nstart = 1,algorithm = c("Lloyd"), trace=FALSE)
#kc3_3 <- kmeans(data_to_plot, 3, iter.max = 50, nstart = 1,algorithm = c("MacQueen"), trace=FALSE)

#kc1_5 <- kmeans(data_to_plot, 5, iter.max = 50, nstart = 1,algorithm = c("Hartigan-Wong"), trace=FALSE)
#kc2_5 <- kmeans(data_to_plot, 5, iter.max = 50, nstart = 1,algorithm = c("Lloyd"), trace=FALSE)
#kc3_5 <- kmeans(data_to_plot, 5, iter.max = 50, nstart = 1,algorithm = c("MacQueen"), trace=FALSE)

library(pheatmap)
k=4
clusters <- kmeans(data_to_plot, k, iter.max = 50)
clusters_samples <- kmeans(t(data_to_plot), 3, iter.max = 50)


print(ncol(data_to_plot))
print(nrow(data_to_plot))
print(clusters)

print("cond")
print(cond.col)
print("names kmneans")
names(clusters$cluster)

print("names data")
names(data_to_plot)

require(RColorBrewer)

col2 <- colorRampPalette(brewer.pal(4,"Set2"));

my_colors <- list( 
          condition = colorRampPalette(c('#FF0000', '#000000', '#0000FF'))(length(cond_continuous))[rank(cond_continuous)]
)

#my_colors <- list( 
#          condition = col2(length(cond_continuous))[rank(cond_continuous)]
#)

names(my_colors$condition) <- rank(cond_continuous)


#my_colors = list(
#    condition = c("Low SSQ" = "#FF0000", "Mid SSQ" = "#000000", "High SSQ" = "#0000FF")
#    
#)



m.kmeans<- cbind(data_to_plot, clusters$cluster) # combine the cluster with the matrix
m.kmeans_samples<- cbind(t(data_to_plot), clusters_samples$cluster) # combine the cluster with the matrix

# the last column
o<- order(m.kmeans[,ncol(m.kmeans)]) # order the last column
o_samples<- order(m.kmeans_samples[,ncol(m.kmeans_samples)]) # order the last column

m.kmeans<- m.kmeans[o,] # order the matrix according to the order of the last column
m.kmeans<- m.kmeans[,o_samples] # order the matrix according to the order of the last column

my_colors$condition <- my_colors$condition[o_samples]
names(my_colors$condition) <- cond_continuous[o_samples]


my_colors <- list( 
          condition = colorRampPalette(c('#FF0000', '#000000', '#0000FF'))(3),
          category = c("Low SSQ" = "#FF0000", "Mid SSQ" = "#000000", "High SSQ" = "#0000FF"),
          raw = colorRampPalette(c('#FF0000', '#000000', '#0000FF'))(3)
          

)

my_colors_subset <- list( 
          condition = colorRampPalette(c('#FF0000', '#0000FF'))(2),
          category = c("Low SSQ" = "#FF0000", "Mid SSQ" = "#000000", "High SSQ" = "#0000FF"),
          raw = colorRampPalette(c('#FF0000', '#0000FF'))(2)
)

#my_sample_col <- data.frame(condition = cond_continuous[o_samples])
my_sample_col <- data.frame(condition = cond_continuous, category=cond, raw=cond_continuous_raw)
my_sample_col <- my_sample_col[o_samples,]

#my_sample_col <- data.frame(condition = cond[o_samples])
row.names(my_sample_col) <- colnames(m.kmeans)


print(str(clusters))
annotation_row = data.frame(feature_cluster = clusters$cluster[o])
#annotation_col_cluster = data.frame(sample_cluster = clusters_samples$cluster[o_samples], condition = cond[o_samples])

annotation_col_cluster = data.frame(sample_cluster = clusters_samples$cluster,condition = cond_continuous, category=cond, raw=cond_continuous_raw)

print("SAMPLE COL UNCLUSTERED")
print(annotation_col_cluster)


annotation_col_cluster = annotation_col_cluster[o_samples,]


print("SAMPLE COL CLUSTERED")
print(annotation_col_cluster)


subset <- which(cond[o_samples]!="Mid SSQ")
m.kmeans_subset <- m.kmeans[,subset]
my_sample_col_subset <- my_sample_col[subset,]
annotation_col_cluster_subset <- annotation_col_cluster[subset,]

str(my_colors)
str(annotation_col_cluster)
str(my_sample_col)

fname=paste0("fig4_normalized_Zscores_cliped_correl_kmeans_dendo_", suffix, ".pdf")


#p <- pheatmap(data_to_plot[names(clusters$cluster), ],
p <- pheatmap(m.kmeans[,1:(ncol(m.kmeans)-1)],
        cluster_rows = FALSE, cluster_cols = TRUE, 
        color=colpal, annotation_col = my_sample_col, annotation_colors = my_colors,
        #color=colpal, annotation_col = my_sample_col,
        filename = fname, height=6,width=12, border_color='NA',
        annotation_row = annotation_row
     )

fname=paste0("fig4_normalized_Zscores_cliped_correl_kmeans_kmeans_", suffix, ".pdf")


#p <- pheatmap(data_to_plot[names(clusters$cluster), ],
p <- pheatmap(m.kmeans[,1:(ncol(m.kmeans)-1)],
        cluster_rows = FALSE, cluster_cols = FALSE, 
        color=colpal, annotation_col = annotation_col_cluster, annotation_colors = my_colors,
        #color=colpal, annotation_col = annotation_col_cluster,
        filename = fname, height=6,width=12, border_color='NA',
        annotation_row = annotation_row
     )

fname=paste0("fig4_normalized_Zscores_cliped_correl_kmeans_kmeans_subset_", suffix, ".pdf")


#p <- pheatmap(data_to_plot[names(clusters$cluster), ],
p <- pheatmap(m.kmeans_subset[,1:(ncol(m.kmeans_subset)-1)],
        cluster_rows = FALSE, cluster_cols = FALSE, 
        color=colpal, annotation_col = annotation_col_cluster_subset, annotation_colors = my_colors_subset,
        #color=colpal, annotation_col = annotation_col_cluster,
        filename = fname, height=6,width=12, border_color='NA',
        annotation_row = annotation_row
     )

fname=paste0("kmeans_centroids_dendo_", suffix, ".pdf")

p <- pheatmap(m.kmeans[,1:(ncol(m.kmeans)-1)],
        kmeans_k = 4, 
        color=colpal, annotation_col = my_sample_col, annotation_colors = my_colors,
        filename = fname, height=6,width=12, border_color='NA',
     )


fname=paste0("fig4_normalized_Zscores_cliped_correl_clustering_complete_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("correlation", "correlation"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE, clust.method = c("complete","complete")
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()

fname=paste0("fig4_normalized_Zscores_cliped_correl_clustering_single_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("correlation", "correlation"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE, clust.method = c("single","single")
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()

fname=paste0("fig4_normalized_Zscores_cliped_correl_clustering_median_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("correlation", "correlation"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE, clust.method = c("median","median")
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()
fname=paste0("fig4_normalized_Zscores_cliped_correl_clustering_average_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("correlation", "correlation"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE, clust.method = c("average","average")
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()
fname=paste0("fig4_normalized_Zscores_cliped_correl_clustering_mcquitty_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("correlation", "correlation"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE, clust.method = c("mcquitty","mcquitty")
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()

fname=paste0("fig4_normalized_Zscores_cliped_correl_clustering_centroid_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("correlation", "correlation"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE, clust.method = c("centroid","centroid")
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()

####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################

#the agglomeration method to be used. This should be (an unambiguous abbreviation of) one of 
#"ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).

fname=paste0("fig4_normalized_Zscores_cliped_euclidean_clustering_ward_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE, clust.method = c("ward","ward")
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()

fname=paste0("fig4_normalized_Zscores_cliped_euclidean_clustering_complete_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE, clust.method = c("complete","complete")
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()

fname=paste0("fig4_normalized_Zscores_cliped_euclidean_clustering_single_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE, clust.method = c("single","single")
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()

fname=paste0("fig4_normalized_Zscores_cliped_euclidean_clustering_median_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE, clust.method = c("median","median")
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()
fname=paste0("fig4_normalized_Zscores_cliped_euclidean_clustering_average_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE, clust.method = c("average","average")
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()
fname=paste0("fig4_normalized_Zscores_cliped_euclidean_clustering_mcquitty_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE, clust.method = c("mcquitty","mcquitty")
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()

fname=paste0("fig4_normalized_Zscores_cliped_euclidean_clustering_centroid_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=names(dataset[,start:ncol(controls)]),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE, clust.method = c("centroid","centroid")
   #lhei=c(1,0.1,1), scale=TRUE
    )
dev.off()


##############################################################
library("RColorBrewer")

library("gplots")

col <- colpal
#col<- colorRampPalette(c("red", "white", "blue"))(256)

#pdf(file="fig4_normalized_z-scores_RD_heatmaps.pdf")

#heatmap(data_to_plot, scale = "none", col =  col, 
#        ColSideColors = cond.col[cond] )

#heatmap.2(scale(data_to_plot), scale = "none", col =  col, 
#        ColSideColors = cond.col[cond], 
#        trace = "none", density.info = "none",
#        margins = c(5, 10)
#        )

#dev.off()
#        RowSideColors = rep(c("blue", "pink"), each = 16),



pca_data <- as.matrix(dataset[,start:ncol(dataset)])

cond=as.character(dataset[,2])
cond.col <- c("Low SSQ" = "red", "Mid SSQ" = "black", "High SSQ" = "blue")


my_pca <- prcomp(pca_data, center = TRUE)

library(ggbiplot)

ggbiplot(my_pca,ellipse=TRUE,  labels=dataset[,1], groups=cond)

library("factoextra")

eig.val <- get_eigenvalue(my_pca)


fviz_eig(my_pca, addlabels = TRUE, ylim = c(0, 20))

library("corrplot")
var <- get_pca_var(my_pca)
corrplot(var$cos2[,1:20], is.corr=FALSE, tl.col = "black")
corrplot(var$contrib[,1:20], is.corr=FALSE, tl.col = "black")
corrplot(var$cor[,1:20], is.corr=TRUE, tl.col = "black")
str(var)
var$cos2

fviz_pca_ind(my_pca, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = cond, 
             col.ind = "black", 
             palette = "jco", 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Condition")


fviz_pca_var(my_pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )

eig.val
ndim=20
data_to_plot <- my_pca$x[,1:ndim]


for(i in c(1:ndim) ){

    a <- 0
#<- scale(dataset[,i],center=TRUE,scale=TRUE)
    data_to_plot[,i] <- scale(data_to_plot[,i],center=TRUE,scale=TRUE)
    #dataset[,i] <- boxCoxVariable(dataset[,i])
    #print(min(dataset[,i]))

}

for(i in 1:nrow(data_to_plot)){

    data_to_plot[i,] <- myClip( data_to_plot[i,] )

}


data_to_plot = t(as.matrix(data_to_plot))


str(my_pca)

my_pca$x

nrow(my_pca$x)

PCA <- princomp(pca_data,cor=T)
PCA
PCA$loadings

fname=paste0("fig5_PCA_", suffix, ".pdf")

pdf(file=fname, height=6,width=12)

cim(data_to_plot, col.sideColors = cond.col[cond], col.names=dataset[,1], row.names=c(1:ndim),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "euclidean"),
    #legend=list( legend = names(cond.col), col = cond.col), dist.method = c("euclidean", "correlation"),
    legend=list( legend = names(cond.col), col = cond.col), dist.method = c("correlation", "correlation"),
    margins = c(5, 13), keysize=c(1.25,1), color=colpal, scale=TRUE
   #lhei=c(1,0.1,1), scale=TRUE
    )

dev.off()


library(FactoMineR)
my_pca <- PCA(pca_data, ncp = 3, graph = FALSE)

res.hcpc <- HCPC(my_pca, graph = FALSE)

fviz_dend(res.hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.8      # Augment the room for labels
          )

fviz_cluster(res.hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
             )

q()

library(Rtsne)

tsne <- Rtsne(pca_data, dims = 2, perplexity=10, verbose=TRUE, max_iter = 50000)

plot(tsne$Y, t='n', main="tsne")
text(tsne$Y, labels=dataset[,1], col=cond.col)

tsne <- Rtsne(pca_data, check_duplicates = FALSE, pca = FALSE, perplexity=10, theta=0.5, dims=2,  max_iter = 50000)

#display results of t-SNE
require(rgl)
plot3d(tsne$Y, col=cond.col)
legend3d("topright", legend = '0':'5', pch = 16, col = cond.col)

rgl.postscript("t-sne_persp3dd.pdf","pdf") 





library(ggplot2)
tsne_plot <- data.frame(x = tsne$Y[,1], y = tsne$Y[,2], col = cond)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))


