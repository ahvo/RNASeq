# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 12:54:02 2016

@author: Andy
"""

import os
#Might need rename files first: for f in *_1.fastq; do mv $f $(echo ${f} |sed  's/_1.fastq/_R1_.fastq.gz/'); done
InputDir = raw_input("Where are your Fastqs located? (include first but not final slash - ie: '/test/2/Files'):")

#Make Analysis Directories
RootFolder = InputDir+'/RNA_Analysis'
ScriptFolder = InputDir+'/RNA_Analysis/Scripts'
BamFolder = InputDir+'/RNA_Analysis/BAM'
Analysis = InputDir+'/RNA_Analysis/Analysis'
FileList = []
CSVList = 'Comparisons'

if not os.path.exists(RootFolder):
    os.makedirs(RootFolder)
    
if not os.path.exists(ScriptFolder):
    os.makedirs(ScriptFolder)
    
if not os.path.exists(BamFolder):
    os.makedirs(BamFolder)
    
if not os.path.exists(Analysis):
    os.makedirs(Analysis)

listdir = os.listdir(InputDir)

#Make Downstream Analysis Scripts (Combine htseqcounts and run them through edgeR) Need to run this script AFTER alignment/counts.  Just run yourself without super server
analysis_script = open(InputDir+'/RNA_Analysis/Scripts/analysis.txt', "w")
print >> analysis_script, 'module load R/3.3.1'
print >> analysis_script, 'perl /projects/p20742/tools/makeHTseqCountsTable.pl', BamFolder, '/projects/p20742/anno/Ens/mm10.Ens_78/mm10.Ens_78.cuff.gtf', Analysis
print >> analysis_script, 'Rscript /projects/p20742/tools/runEdgeRrnaSeq.2.R --assembly=mm10 --countFile='+Analysis+'/htseq.all.counts.txt', '--numCores=4 --outputDirectory='+Analysis, '--runMDS=1'
print >> analysis_script, 'Rscript /projects/p20742/tools/runEdgeRrnaSeq.2.R --assembly=mm10 --countFile='+Analysis+'/htseq.all.counts.txt', '--comparisonFile='+Analysis+'/comparison.csv --numCores=4 --outputDirectory='+Analysis, '--runMDS=0'
print >> analysis_script, 'for f in', Analysis+'/*.edgeR.txt; do Rscript /projects/p20742/tools/makeMAplot.R --degFile=${f}', '--adjp=0.01 --labelTop=1; done'
print >> analysis_script, 'for f in', Analysis+'/*.edgeR.txt; do Rscript /projects/p20742/tools/makeBigHeatmap.R --degFile=${f}', '--adjp=0.01 --countFile='+Analysis+'/htseq.normCounts.txt; done'
print >> analysis_script, 'for f in', Analysis+'/*.edgeR.txt; do Rscript /projects/p20742/tools/runGOforDEG.R --degFile=${f}', '--adjp=0.01 --assembly=mm10; done'
analysis_script.close()
for file in listdir:
    if file.endswith("R1_.fastq.gz"):
        #Make gene list to comparison CSV
        FileList.append(file[0:len(file)-13])
        #Makes Shell script with file name
        script = open(InputDir+'/RNA_Analysis/Scripts/'+file[0:len(file)-13]+'.sh', "w")
        #Print MSUB Header
        script.write('''#!/bin/bash
#MSUB -A b1042
#MSUB -q genomics
#MSUB -l walltime=24:00:00
#MSUB -m a
#MSUB -j oe
#MSUB -l nodes=1:ppn=4
export PATH=$PATH:/projects/p20742/tools/
module load bowtie2/2.2.6
module load tophat/2.1.0
module load samtools
module load boost
module load bowtie
module load gcc/4.8.3
module load java
module load R/3.3.1

export PATH=$PATH:/projects/b1042/McNallyLab/AV_RNA_ANALYSIS/TOOLS/RSEM-1.3.0/
export PATH=$PATH:/projects/p20742/tools/

''')      
        R1Path = InputDir+'/'+file
        R2Path = InputDir+'/'+file[0:len(file)-11]+'2_.fastq.gz'
        # Run Tophat to align data for paired end data.
        print >> script, 'tophat --no-novel-juncs --read-mismatches 2 --read-edit-dist 2 --num-threads 4 --max-multihits 5 --transcriptome-index /projects/p20742/anno/tophat_tx/mm10.Ens_78.cuff -o', BamFolder+'/'+file[0:len(file)-13], '/projects/p20742/anno/bowtie_indexes/mm10', R1Path, R2Path, '>&', BamFolder+'/'+file[0:len(file)-13]+'.tophat.log'
        print >> script, 'mv', BamFolder+'/'+file[0:len(file)-13]+'/accepted_hits.bam', BamFolder+'/'+file[0:len(file)-13]+'.bam'   
        print >> script, 'samtools sort -n', BamFolder+'/'+file[0:len(file)-13]+'.bam', BamFolder+'/'+file[0:len(file)-13]+'_NameSorted'
        print >> script, 'module unload mpi'
        print >> script, 'module load python'
        print >> script, 'htseq-count -f bam -q -m intersection-nonempty -s no -t exon -i gene_id', BamFolder+'/'+file[0:len(file)-13]+'_NameSorted.bam', '/projects/p20742/anno/Ens/mm10.Ens_78/mm10.Ens_78.cuff.gtf >', BamFolder+'/'+file[0:len(file)-13]+'.htseq.counts'
        print >> script, 'module unload python' 
        print >> script, 'module load gcc/4.8.3'
        #print >> script, 'perl /projects/p20742//tools/makeHTseqCountsTable.pl', /projects/b1042/McNallyLab/AV_RNA_ANALYSIS/LV/LV/bam// /projects/p20742//anno/Ens/mm10.Ens_78/mm10.Ens_78.cuff.gtf /projects/b1042/McNallyLab/AV_RNA_ANALYSIS/LV/LV/bam/

        #Apply RSEM
        #print >> script, '/projects/b1042/McNallyLab/AV_RNA_ANALYSIS/TOOLS/RSEM-1.3.0/rsem-calculate-expression --paired-end', InputDir+'/'+UnzipName, InputDir+'/'+file[0:len(file)-11]+'3_.fastq', '/projects/p20742/anno/rsemTx/mm10.Ens_78', InputDir+'/RNA_Analysis/BAM/'+file[0:len(file)-13], '--no-bam-output -p 4  --estimate-rspd >&', InputDir+'/RNA_Analysis/BAM/'+file[0:len(file)-13]+'.rsem.log'+'\n'
        #print >> script, '/projects/b1042/McNallyLab/AV_RNA_ANALYSIS/TOOLS/RSEM-1.3.0/rsem-plot-model', InputDir+'/RNA_Analysis/BAM/'+file[0:len(file)-13], InputDir+'/RNA_Analysis/BAM/'+file[0:len(file)-13]+'.rsemPlot.pdf'+'\n'
        #Make RSEM counts Table
        #print >> script, 'perl /projects/p20742/tools/makeRSEMcountsTable.pl', InputDir+'/RNA_Analysis/BAM/', '/projects/p20742/anno/Ens/mm10.Ens_78/mm10.Ens_78.cuff.gtf', InputDir+'/RNA_Analysis/BAM/', 'isoforms'+'\n'
        #print >> script, 'Rscript /projects/p20742/tools/runEdgeRrnaSeq.2.R --assembly=mm10 --countFile='+InputDir+'/RNA_Analysis/BAM/rsem.all.counts.txt --numCores=4 --outputDirectory='+InputDir+'/RNA_Analysis --runMDS=1'+'\n'
        #Close Script File
        script.close()

#make comparison File (1 is reference and -1 is experimental)
ComparisonCSV = open(Analysis+'/comparison.csv', "w")    
for title in FileList:
    CSVList = CSVList+','+ title

print >> ComparisonCSV, CSVList+'\n'+'Input Comparison Here,1,0,-1'+'\n'

ComparisonCSV.close()
 
        
        
        
        
        
        
