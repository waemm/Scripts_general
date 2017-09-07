import os
import sys
import subprocess


raw_read_dir = "/SAN/biomed/biomed1/ingest2" #sys.argv[1]

out_dir = sys.argv[2]

cpu_ = sys.argv[3]

name = sys.argv[4]

read_ending = sys.argv[5]

ref = "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/STAR" #sys.argv[6]

id_sep = "_R" #sys.argv[7]

id_index = 0 #int(sys.argv[8])

#mate_pref = sys.argv[7]

#mate_ind = sys.argv[8]

#mate_sep = 

###########################################


star_engine = "/SAN/biomed/biomed6/biomed6/warren/Applications/STAR/update_STAR/STAR/bin/Linux_x86_64_static/STAR"

if ref == "hg":
    print "USING HG19_v2.3 !"
    ref ="/SAN/biomed/biomed10/ingest1/warren/reference/STAR_2.4_updated_hg19_gencode19" # "/SAN/biomed/biomed6/biomed6/warren/Applications/Star_reference/hg19_v2.3"
#ref = "/SAN/biomed/biomed6/biomed6/warren/Applications/Star_reference/hpv_16"
###########################################


"""
BSSE_QGF_25350_141031_D00535_0030_BC57VRANXX_8_TTAGGC_L008_R 1_001_Ctl_a.fastq.gz
BSSE_QGF_25350_141031_D00535_0030_BC57VRANXX_8_TTAGGC_L008_R 1_002_Ctl_a.fastq.gz
BSSE_QGF_25350_141031_D00535_0030_BC57VRANXX_8_TTAGGC_L008_R 2_001_Ctl_a.fastq.gz
BSSE_QGF_25350_141031_D00535_0030_BC57VRANXX_8_TTAGGC_L008_R 2_002_Ctl_a.fastq.gz
"""


read_pairs = {}

outputlist = []
for (path, dirs, files) in os.walk(raw_read_dir):
    
    print (path, dirs, files)
    for readfile in files:
        
        if readfile.endswith(read_ending):
            
            ##### add in sort for paired end files!!!!!!!!!!!!!
            
            root = os.path.split(path)[1].replace(" ","_")

            print "sample type :",root 
            
            #root = id_sep.join(readfile.split(id_sep)[:id_index])
            
            print "root",root
            
            readfile= os.path.join(path,readfile)
            
            #readfile_root = readfile.split(id_sep)[id_index]
            
            readfile_root = readfile.replace("_R1_","_")
            readfile_root = readfile_root.replace("_R2_","_")
            
            if root not in read_pairs:
                
                read_pairs[root] = {}
            
            if readfile_root not in read_pairs[root]: 
                read_pairs[root][readfile_root] = []
                
            read_pairs[root][readfile_root].append(readfile.replace(" ","\ "))

for r in read_pairs:
    
    if len(read_pairs[r]) >= 1:
        
        print r
        
        print len(read_pairs[r])
    
print read_pairs    



    

tmp_str = '#$ -l tmem=5G\n#$ -S /bin/bash\n#$ -l h_vmem=5G\n#$ -l h_rt=98:0:0\n#$ -j y\n#$ -R y\n#$ -pe smp 9\n#$ -cwd\nPATH=/usr/java/latest/jre/bin/:/home/regmwem/local/R/bin:/home/regmwem/local/bin/python2.7/bin:/SAN/biomed/biomed10/ingest1/warren/bif_tools/annovar:/SAN/biomed/biomed6/biomed6/warren/Applications/bowtie2-2.0.5/:/SAN/biomed/biomed6/biomed6/warren/Applications/tophat-2.0.6.Linux_x86_64:/share/apps/genomics/samtools-0.1.18/:$PATH\nif [ ! -d "/scratch0/regmwem" ]; then\n\tmkdir /scratch0/regmwem\nfi\n'
tmp_str+= 'cd %s\n'%(out_dir)
for readp in read_pairs:
    
    first_mates = []
    second_mates = []
    
    if len(read_pairs[readp]) >= 1:
        
        
        print read_pairs[readp]
        
        
        for sample_part in read_pairs[readp]:
            
            samples_ = read_pairs[readp][sample_part]
            
            samples_.sort()
            
            print "order should be 1 2 ",samples_
            
            raw_input()
            first_mates.append(samples_[0])
            second_mates.append(samples_[1])
        
        
        
        read1 = ",".join(first_mates)

        read2 = ",".join(second_mates)
        
        log = "%s.sort.bam"%(readp)
        outname = "%s_out"%(readp)
        
        output_dir_sample = os.path.join(out_dir,"%s_%s"%(readp,name))
        
        
        #### --readFilesIn Read1a.gz,Read1b.gz Read2a.gz,Read2b.gz
        
        #tmp_str +="mkdir %s\ncd %s\n"%(output_dir_sample,output_dir_sample)
        if len(read_pairs[readp]) > 1:
            tmp_str +="%s --genomeDir %s --twopassMode Basic --outSAMunmapped Within --chimSegmentMin 15 --readFilesIn %s %s --readFilesCommand zcat --limitSjdbInsertNsj 2000000 --limitOutSJcollapsed 2000000 --runThreadN %s --genomeLoad NoSharedMemory  --outFileNamePrefix %s\n"%(star_engine,ref,read1,read2,cpu_,"%s_%s"%(readp,name))
        else:
            tmp_str +="%s --genomeDir %s --twopassMode Basic --outSAMunmapped Within --chimSegmentMin 15 --readFilesIn %s --readFilesCommand zcat --runThreadN %s --limitSjdbInsertNsj 2000000 --limitOutSJcollapsed 2000000 --genomeLoad NoSharedMemory --outFileNamePrefix %s\n"%(star_engine,ref,read1,cpu_,"%s_%s"%(readp,name))
        
tmp_str +="cd ..\n%s --genomeLoad Remove"%(star_engine)
        
open('run_star_all_%s.sh'%(name),'w').write(tmp_str)
