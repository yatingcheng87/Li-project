args = c("pten","pten.csv")

# Get counts file from analysis/fname/fname.T.csv
bname = basename(args[1]) 
fname = paste(bname,"-count.T.csv",sep='')
dat = read.csv(file.path("analysis",bname,fname), header = TRUE, row.names=1)
keep = rowSums(cpm(dat) > 3) >= 3
counts = dat[keep,]


#sel = grepl("MT-.*", rownames(counts)) + grepl("ERCC-.*", rownames(counts)) + grepl("mt-.*", rownames(counts))

#counts = counts[!sel,]


head(counts)
factors <- c("dko", "sko", "sko", "sko", "dko", "dko","dko" ,"sko")
groups <- factors

y = DGEList(counts=counts, group=factors)
y = calcNormFactors(y)
f=y$samples$lib.size * y$samples$norm.factors/mean(y$samples$lib.size*y$samples$norm.factors)
y = estimateCommonDisp(y)
y = estimateTagwiseDisp(y)



scaled.counts=round(t(t(counts)/f)*mean(f))
rownames(scaled.counts) = rownames(counts)
dfs = split.data.frame(t(scaled.counts), groups)
dfss = sapply(dfs, colMeans)
lrt = exactTest(y, pair = c(2,1))
ot1 = topTags(lrt,n=nrow(counts),sort.by="PValue")$table
topTags(lrt, n =10)


ot1 = merge(ot1, dfss, by=0)
ot1 = ot1[order(ot1$FDR),] # Sort by ascending FDR
write.csv(ot1,"dko_vssko-yating.csv",row.names=FALSE)
p = plotMDS(y)
p = plotBCV(y)
