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

#id_index = int(sys.argv[8])
#print map(int,sys.argv[8].split(","))
id_index = map(int,sys.argv[8].split(","))
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



read_pairs = {}

outputlist = []
for (path, dirs, files) in os.walk(raw_read_dir):
    
    for readfile in files:
        
        if readfile.endswith(read_ending):
            
            ##### add in sort for paired end files!!!!!!!!!!!!!
            
            #root = os.path.split(path)[1]

            print readfile 
	    #print readfile.split(id_sep)[id_index]
            root = "_".join(readfile.split(id_sep)[id_index[0]:id_index[1]])
            
            print "root",root
            
            
            readfile= os.path.join(path,readfile)
            
            if root not in read_pairs:
                
                read_pairs[root] = [readfile]
                
            else:
                
                read_pairs[root].append(readfile)

for r in read_pairs:
    
    if len(read_pairs[r]) >= 1:
        
        print r
    
print read_pairs    
    

tmp_str = '#$ -l tmem=5G\n#$ -S /bin/bash\n#$ -l h_vmem=5G\n#$ -l h_rt=98:0:0\n#$ -j y\n#$ -R y\n#$ -pe smp 9\n#$ -cwd\nPATH=/usr/java/latest/jre/bin/:/home/regmwem/local/R/bin:/home/regmwem/local/bin/python2.7/bin:/SAN/biomed/biomed10/ingest1/warren/bif_tools/annovar:/SAN/biomed/biomed6/biomed6/warren/Applications/bowtie2-2.0.5/:/SAN/biomed/biomed6/biomed6/warren/Applications/tophat-2.0.6.Linux_x86_64:/share/apps/genomics/samtools-0.1.18/:$PATH\nif [ ! -d "/scratch0/regmwem" ]; then\n\tmkdir /scratch0/regmwem\nfi\n'
tmp_str+= 'cd %s\n'%(out_dir)
for readp in read_pairs:
    
    if len(read_pairs[readp]) >= 1:
        
        read1 = os.path.split(read_pairs[readp][0])[1]
	
	
	
	
	
        #tmp1 = os.path.join(raw_read_dir,read_pairs[readp][0])+'.sai'
        if len(read_pairs[readp]) > 1:
            read2 = os.path.split(read_pairs[readp][1])[1]
        #tmp2 = os.path.join(raw_read_dir,read_pairs[readp][1])+'.sai'
        
        #bam1 = os.path.join(output_dir,readp+'.paired.bam')
        #os.path.join(output_dir,readp+'.paired_sorted')
        
        log = "%s.sort.bam"%(readp)
        outname = "%s_out"%(readp)
        
        output_dir_sample = os.path.join(out_dir,"%s_%s"%(readp,name))
        
        # /SAN/biomed/biomed6/biomed6/warren/Applications/trim_galore_zip/trim_galore
        #wt-HAP-ribotag-1_trimmed.fq.gz
        # --length 15
        
        
        #tmp_str +="mkdir %s\ncd %s\n"%(output_dir_sample,output_dir_sample)
        if len(read_pairs[readp]) > 1:
	    tmp_str += "/SAN/biomed/biomed6/biomed6/warren/Applications/TrimGalore/trim_galore --length 14 --paired %s %s\n"%(read1,read2)
            tmp_str +="%s --genomeDir %s --twopassMode Basic  --readFilesIn %s_trimmed.fq.gz %s_trimmed.fq.gz --readFilesCommand zcat --limitSjdbInsertNsj 2000000 --limitOutSJcollapsed 2000000 --runThreadN %s --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --outFileNamePrefix %s\n"%(star_engine,ref,read1.split(".")[0],read2.split(".")[0],cpu_,"%s_%s"%(readp,name))
        else:
	    
	    tmp_str += "/SAN/biomed/biomed6/biomed6/warren/Applications/TrimGalore/trim_galore --length 14 %s\n"%(read1)
	    
            tmp_str +="%s --genomeDir %s --twopassMode Basic  --readFilesIn %s_trimmed.fq.gz --readFilesCommand zcat --runThreadN %s --limitSjdbInsertNsj 2000000 --limitOutSJcollapsed 2000000 --genomeLoad NoSharedMemory --outSAMtype BAM SortedByCoordinate --outFileNamePrefix %s\n"%(star_engine,ref,read1.split(".")[0],cpu_,"%s_%s"%(readp,name))
        
tmp_str +="cd ..\n%s --genomeLoad Remove"%(star_engine)
        
open('run_star_all_%s.sh'%(name),'w').write(tmp_str)
