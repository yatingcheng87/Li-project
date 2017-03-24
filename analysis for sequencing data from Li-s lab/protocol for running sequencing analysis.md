#### Please check the link for detailed protocol for general analysis procedure #####

https://github.com/chapkinlab/sequencing-pipeline/tree/master/docs

#### Specific procedure running for sequencing data from Dr. Li######

####  raw data downloaded from the website: https://download.txgen.tamu.edu/qinglei.li/160517_D00572_0193_AC8KJKANXX_15375Li/
    
wget https://download.txgen.tamu.edu/qinglei.li/160517_D00572_0193_AC8KJKANXX_15375Li/001/fastq_15375Li_N16029_L001.tar

#### md5sum check 
[cyating@nfsc-oracle-tamu-edu ~]$ md5sum fastq_15375Li_N16029_L001.tar
398d0e0ee403ddb55d007652e6a94dd3  fastq_15375Li_N16029_L001.tar

####  extract files and remove unnessary files in the folder
tar -xvf fastq_15375Li_N16029_L001.tar
rm Undetermined*
rm *.ADAPTER*
rm *.FILTERED*

#### no need to concatenate the files, therefore we can directly gunzip and bzip2 files. 
gunzip -c file.gz | bzip2 -c > file.bz2

ls ~/15375Li_N16029/*.bz2 > listname.txt  # the listname.txt can be used for list file (lists/pten)

#### map reads to the genome
main-scripts/map.sh lists/pten 2> err.log | tee out.log

#### analysis summarization
main-scripts/summary.py lists/pten keys/pten.csv #  you need to create key file before running the summary.py

#### differential analysis will be carried out in the R interface (use edger-yating.R), for convenience, you can transfer the analysis folder, list file, key file to your local computer for running R.

Problem with the current R script
"It looks like two or more levels of classification are needed to run the differential expression script. It seems to work if I add one more test column to your key file. Thanks for finding this. If you're interested, use this key file and add "test" to the factors and design commands, and try. It should work. I'll talk to Jason and fix this. " - from Kumaran. 

"the line staring with scaled.counts in your code is supposed to save the normalized counts. However, the equation for obtaining the normalized counts turned out to be wrong and the line should be changed as follows:
y=calcNormFactors(y)
f=y$samples$lib.size * y$samples$norm.factors/mean(y$samples$lib.size*y$samples$norm.factors)
scaled.counts=round(t(t(counts)/f)*mean(f))" - from Eunji

"sel = grepl("MT-.*", rownames(counts)) + grepl("ERCC-.*", rownames(counts)) + grepl("mt-.*", rownames(counts))
Originally, this line was written in order to remove ERCCs and mitochondria under the assumption that their names start with ERCC and MT, respectively. However, you can easily check out that this code removes other genes as well. Also, if you want to keep everything, you should change this line accordingly or delete the part." - from Eunji

Therefore, if you only have one level of classification, you have to add an additional artificial one. This step might be confusing. 
In addition, you can read the output results when using the R script in the R console. It is easier for you to catch up error and understand each step. 

There are different ways that you can utilize the result you get from the pipeline  in the R. For the detail, please check the detailed comment in the R script. exact-test-for Li-3-23-2017.R. I have incoporated what Eunji suggested in the R script. 
(https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf  The edgeR manual will be good sourse to get help. )

Since there are only two groups for comparision, we used the exact T test, if you have multiple groups or multiple classification, you can easily modify the script to run the GLM test.



If you don't have the edgeR package in your local computer. you need to download the edgeR package

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)


 
