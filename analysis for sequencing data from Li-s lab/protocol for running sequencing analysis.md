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


 