import os
import sys
import subprocess


raw_read_dir =  "/cluster/project9/MML2/FASTQ/" #sys.argv[1]

out_dir = "/cluster/project9/MML2/" # sys.argv[2]

cpu_ = 10 #sys.argv[3]

name = "MML2"#sys.argv[4]

read_ending = "fastq.gz"#sys.argv[5]

ref = "/SAN/vyplab/HuRNASeq/reference_datasets/STAR/mouse/" #sys.argv[6]

id_sep = "."#sys.argv[7]

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



read_pairs = {}

outputlist = []
for (path, dirs, files) in os.walk(raw_read_dir):
    
    for readfile in files:
        
        if readfile.endswith(read_ending):
            
            ##### add in sort for paired end files!!!!!!!!!!!!!
            
            #root = os.path.split(path)[1]

            #print readfile 
	    #print readfile.split(id_sep)[id_index]
            root = os.path.split(path)[1]
            
            print "root",root
            
            
            readfile_fullpath = os.path.join(path,readfile)
            
            if root not in read_pairs:
                
                read_pairs[root] = {"R1":[],"R2":[]}
                
           
		
	    if readfile.find("R1") > -1:
		read_pairs[root]["R1"].append(readfile_fullpath)
		
	    elif readfile.find("R2") > -1:
		read_pairs[root]["R2"].append(readfile_fullpath)
		
		
	    
                
	    #read_pairs[root].append(readfile)


    

tmp_str = '#$ -l tmem=5G\n#$ -S /bin/bash\n#$ -l h_vmem=5G\n#$ -l h_rt=98:0:0\n#$ -j y\n#$ -R y\n#$ -pe smp 9\n#$ -cwd\nPATH=/usr/java/latest/jre/bin/:/home/regmwem/local/R/bin:/home/regmwem/local/bin/python2.7/bin:/SAN/biomed/biomed10/ingest1/warren/bif_tools/annovar:/SAN/biomed/biomed6/biomed6/warren/Applications/bowtie2-2.0.5/:/SAN/biomed/biomed6/biomed6/warren/Applications/tophat-2.0.6.Linux_x86_64:/share/apps/genomics/samtools-0.1.18/:$PATH\nif [ ! -d "/scratch0/regmwem" ]; then\n\tmkdir /scratch0/regmwem\nfi\n'
tmp_str+= 'cd %s\n'%(out_dir)
for readp in read_pairs:
    
    print "readp : ",readp
    
    if len(read_pairs[readp]) >= 1:
	
	
	read1 = read_pairs[readp]["R1"]
	
	read1.sort()
	
	print len(read1)
	print read1
	
	read2 = read_pairs[readp]["R2"]
		
	read2.sort()
	
	print len(read2)
	print read2	
        
        
        #raw_input()
        
        log = "%s.sort.bam"%(readp)
        outname = "%s_out"%(readp)
        
        output_dir_sample = os.path.join(out_dir,"%s"%(readp))
        
        
        
        
        #tmp_str +="mkdir %s\ncd %s\n"%(output_dir_sample,output_dir_sample)
        if len(read_pairs[readp]) > 1:
            tmp_str +="%s --genomeDir %s --twopassMode Basic  --readFilesIn %s %s --readFilesCommand zcat --limitSjdbInsertNsj 2000000 --limitOutSJcollapsed 2000000 --runThreadN %s --genomeLoad NoSharedMemory --outSAMtype BAM Unsorted --outFileNamePrefix %s\n"%(star_engine,ref,",".join(read1),",".join(read2),cpu_,"%s_%s"%(readp,name))
       
tmp_str +="cd ..\n%s --genomeLoad Remove"%(star_engine)
        
open('run_star_all_%s.sh'%(name),'w').write(tmp_str)
