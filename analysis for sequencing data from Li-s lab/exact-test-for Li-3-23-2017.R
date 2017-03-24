
###input count data and group annoation in to the R. There are many other ways to import 
###data. In this script, we tried to use the most of command from the original script in the 
###pipeline

## for instance, you can just use the pten-count.T.csv alone to extract count data##

#dat1 = read.csv(file = "pten-count.T.csv", header = TRUE)
#names=dat1[,1]
#rownames(dat1)=make.names(names,unique=TRUE)
#dat <- dat1[,-1] #this deletes the first column 
---------------------------------------------------------------------------------------------------------
### or you can download the analysis file, the list file, key file from the pipeline and use 
###  the original script in the pipeline 
  
args = c("pten","pten.csv")

# Get counts file from analysis/fname/fname.T.csv
bname = basename(args[1]) 
fname = paste(bname,"-count.T.csv",sep='')
dat = read.csv(file.path("analysis",bname,fname), header = TRUE, row.names=1)

#### filter the data ###
keep = rowSums(cpm(dat) > 3) >= 3
counts = dat[keep,]

#### the following two commonds are used to remove ERCC count and mitochondrial count
### use with caution, since it might remove your interested genes
### for instance, the mt- genes might not be mitochodrial genes. 

#sel = grepl("MT-.*", rownames(counts)) + grepl("ERCC-.*", rownames(counts)) + grepl("mt-.*", rownames(counts))

#counts = counts[!sel,]


head(counts)

# there are different ways to make group for the files. 
#first method

factors <- c("dko", "sko", "sko", "sko", "dko", "dko","dko" ,"sko") # the order should match to the order of sample names
#in the count csv file. 

groups <- factors
# second method, generate the group file, group.csv. Be sure the order of sample list 
#is the same as in the count csv file
#metadata = read.csv("group.csv", header=TRUE, row.names=1)
#factors <- metadata$treatment
#groups <- factors

### edgeR script, you can check the manual of edgeR for more detailed algorithm

y = DGEList(counts=counts, group=factors)
y = calcNormFactors(y)
f=y$samples$lib.size * y$samples$norm.factors/mean(y$samples$lib.size*y$samples$norm.factors)
y = estimateCommonDisp(y)
y = estimateTagwiseDisp(y)


## normalized count calculation
scaled.counts=round(t(t(counts)/f)*mean(f))
rownames(scaled.counts) = rownames(counts)
dfs = split.data.frame(t(scaled.counts), groups)
dfss = sapply(dfs, colMeans)
lrt = exactTest(y, pair = c(2,1)) # the edgeR consider the group in alphabetical order
#therefore, group 1 is dko, group 2 is sko, if you want to compare dko vs sko
# pair = c(2,1) is group 1 vs group2, namely, dko vs sko.

ot1 = topTags(lrt,n=nrow(counts),sort.by="PValue")$table
topTags(lrt, n =10)


ot1 = merge(ot1, dfss, by=0)
ot1 = ot1[order(ot1$FDR),] # Sort by ascending FDR
write.csv(ot1,"dko_vssko-yating-modified.csv",row.names=FALSE)
p = plotMDS(y)
p = plotBCV(y)
