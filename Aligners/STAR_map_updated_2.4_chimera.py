import os
import sys
import subprocess


raw_read_dir = sys.argv[1]

out_dir = sys.argv[2]

cpu_ = sys.argv[3]

name = sys.argv[4]

read_ending = sys.argv[5]

ref = sys.argv[6]

id_sep = sys.argv[7]

id_index = int(sys.argv[8])

#mate_pref = sys.argv[7]

#mate_ind = sys.argv[8]

#mate_sep = 

###########################################


star_engine = "/SAN/biomed/biomed6/biomed6/warren/Applications/STAR/bin/Linux_x86_64_static/STAR"

if ref == "hg":
    print "USING HG19_v2.3 !"
    ref = "/SAN/biomed/biomed6/biomed6/warren/Applications/Star_reference/hg19_v2.3"
#ref = "/SAN/biomed/biomed6/biomed6/warren/Applications/Star_reference/hpv_16"
###########################################



read_pairs = {}

####### Collect all read files that correspond to read_ending and put in a hashtable/dictionary

outputlist = []
for (path, dirs, files) in os.walk(raw_read_dir):
    
    for readfile in files:
        
        if readfile.endswith(read_ending):
            
            ##### add in sort for paired end files!!!!!!!!!!!!!
            
            root = readfile.strip().split(id_sep)[id_index]
            
            readfile= os.path.join(path,readfile)
            
            if root not in read_pairs:
                
                read_pairs[root] = [readfile]
                
            else:
                
                read_pairs[root].append(readfile)

for r in read_pairs:
    
    if len(read_pairs[r]) >= 1:
        
        print r
    
print read_pairs    
    

# script start, I ask for 5G x 12 cores = 60GB total

tmp_str = '#$ -l scr=40G\n#$ -l tmem=5G\n#$ -S /bin/bash\n#$ -l h_vmem=5G\n#$ -l h_rt=98:0:0\n#$ -j y\n#$ -R y\n#$ -pe smp 12\n#$ -cwd\nPATH=/usr/java/latest/jre/bin/:/home/regmwem/local/R/bin:/home/regmwem/local/bin/python2.7/bin:/SAN/biomed/biomed10/ingest1/warren/bif_tools/annovar:/SAN/biomed/biomed6/biomed6/warren/Applications/bowtie2-2.0.5/:/SAN/biomed/biomed6/biomed6/warren/Applications/tophat-2.0.6.Linux_x86_64:/share/apps/genomics/samtools-0.1.18/:$PATH\nif [ ! -d "/scratch0/regmwem" ]; then\n\tmkdir /scratch0/regmwem\nfi\n'

for readp in read_pairs:
    
    if len(read_pairs[readp]) >= 1:
        
        read1 = read_pairs[readp][0]
        #tmp1 = os.path.join(raw_read_dir,read_pairs[readp][0])+'.sai'
        if len(read_pairs[readp]) > 1:
            read2 = read_pairs[readp][1]
        #tmp2 = os.path.join(raw_read_dir,read_pairs[readp][1])+'.sai'
        
        #bam1 = os.path.join(output_dir,readp+'.paired.bam')
        #os.path.join(output_dir,readp+'.paired_sorted')
        
        log = "%s.sort.bam"%(readp)
        outname = "%s_out"%(readp)
        
        output_dir_sample = os.path.join(out_dir,"%s_%s"%(readp,name))
        
        
        ########## for my pipeline i create a new folder for each sample, this is probably not the most efficient, use the --outFileNamePrefix
	###### All known junctions where built into the index
	###### I do chimeric read detection but this can be switched off if not needed 
	
	###### This command: 1. unzips the FASTQ on the fly 2. maps to the genome + chimeric reads 3. Co-ordinate sorts reads once done
        
        tmp_str +="mkdir %s\ncd %s\n"%(output_dir_sample,output_dir_sample)
        if len(read_pairs[readp]) > 1:
            tmp_str +="%s --genomeDir %s --readFilesIn %s %s --chimSegmentMin 20  --readFilesCommand zcat --runThreadN %s --genomeLoad LoadAndKeep --outSAMstrandField intronMotif --limitBAMsortRAM 15000000000 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within\n"%(star_engine,ref,read1,read2,cpu_)
        else:
            tmp_str +="%s --genomeDir %s --readFilesIn %s --readFilesCommand zcat --runThreadN %s --genomeLoad LoadAndKeep --outSAMstrandField intronMotif --limitBAMsortRAM 15000000000 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within 2>&1 > %s\n"%(star_engine,ref,read1,cpu_,log)
        
tmp_str +="cd ..\n%s --genomeLoad Remove"%(star_engine)
        
open('run_star_all.sh','w').write(tmp_str)
