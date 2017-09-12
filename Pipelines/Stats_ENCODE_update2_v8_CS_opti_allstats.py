# filter and sort the junctions for TSC2
import os,sys,csv,pysam,numpy,math
from operator import itemgetter
import pybedtools,numpy
#EOI = [[33460189,33460260] ] #   ,[[2127598,2127727] ]   # 

import scipy.stats as stats
import msgpack
import cPickle as pickle


control_id = pickle.load(open("allRBPS/control_id_v4_allENC.p"))


#### number of reads to consider the sample data
minimum_Reads = 7

data_ = {}

count = 0

data_doubles = {}


### ENCODE data
input_table = "/SAN/biomed/biomed10/ingest1/warren/PhylopScores/annotate_exons/Integrate_tables/ENCODEdata_summary_dictionary.csv" # sys.argv[1]# "Lorea_all12_ENCODE_seq_IDs_Update2_junction_data_skipped_junction_data_v1.csv" #"Lorea_all12_ENCODE_seq_IDs_Update2_K562_junction_data_skipped_junction_data_v1.csv" #"Lorea_all10_ENCODE_seq_IDs_skipped_junction_data_v1.csv" #'skip_enc_short_intervals.csv' # ENCODE data "Lorea_all11_ENCODE_seq_IDs_K562_junction_data_skipped_junction_data_v1.csv" #




maxent_file = "ENCODEv2.1_5p_maxent_in.fa.out"# "TEST_ENCODEv4_junction_data_maxent_in.fa.out" #sys.argv[2] # "Lorea_all12_ENCODE_seq_IDs_ENCODEv2_junction_data_maxent_in.fa.out" #"""Lorea_all12_ENCODE_seq_IDs_ENCODEv1_junction_data_maxent_in.fa.out" #"Lorea_all12_ENCODEv2_seq_IDs_ENCODEv2_junction_data_maxent_in.fa.out"

name =  "ENCODE_skip_V7ranksum" #sys.argv[3]






samples= []

control_sample = ""






######## statistics to include

# 'CodonFrame','diff_score','downstream_intron_length', 'Exon_upstream_cons',  'Upstream_unaligned_ESS ', 'Upstream_aligned_ESS ', 'Upstream_aligned_RescueESS ', 'Downstream_unaligned_ESS ', 'Downstream_aligned_ESS ', 'Downstream_aligned_RescueESS '

#stats_fields = ['Low_complexity_number', 'Low_complexity_ratio', 'LTR_number', 'LTR_ratio', 'Satellite_number', 'Satellite_ratio', 'rRNA_number', 'rRNA_ratio', 'DNA_number', 'DNA_ratio', 'Simple_repeat_number', 'Simple_repeat_ratio', 'Unknown_number', 'Unknown_ratio', 'scRNA_number', 'scRNA_ratio', 'RC_number', 'RC_ratio', 'tRNA_number', 'tRNA_ratio', 'RNA_number', 'RNA_ratio', 'LINE_number', 'LINE_ratio', 'SINE_number', 'SINE_ratio', 'snRNA_number', 'snRNA_ratio', 'exon_length']
#stats_fields = []


#stats_fields = ["MIR_distance",'L1_distance','L2_distance','TcMar-Tigger_distance','Alu_distance','exon_length','LINE_ratio', 'SINE_ratio','LTR_ratio','MIR_ratio','L2_ratio','L1_ratio','TcMar-Tigger_ratio','Alu_ratio','Downstream_unaligned_Pentamer ', 'Upstream_aligned_Pentamer ', 'downstream_intron_length', 'Exon_end_cons', 'Exon_start_cons', 'upstream_intron_length']


stats_fields = ['exon_length','Alu_ratio','LTR_ratio','MIR_ratio','L2_ratio','L1_ratio','TcMar-Tigger_ratio','Rs_score']



ranges_ = {"exon_length":[10,20,30,40,50,60,80,100,150,200,250,300,350,400,500,600],
           "MIR_ratio":[ 0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,  0.45,  0.5 , 0.55,  0.6 ,  0.7,  0.8,  0.9],
            "L2_ratio":[ 0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,  0.45,  0.5 , 0.55,  0.6 ,  0.7,  0.8,  0.9],
             "L1_ratio":[ 0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,  0.45,  0.5 , 0.55,  0.6 ,  0.7,  0.8,  0.9],
           "TcMar-Tigger_ratio":[ 0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,  0.45,  0.5 , 0.55,  0.6 ,  0.7,  0.8,  0.9],
           "Alu_ratio":[ 0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,  0.45,  0.5 , 0.55,  0.6 ,  0.7,  0.8,  0.9],
           "LINE_ratio":[ 0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,  0.45,  0.5 , 0.55,  0.6 ,  0.7,  0.8,  0.9],
           "SINE_ratio":[ 0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,  0.45,  0.5 , 0.55,  0.6 ,  0.7,  0.8,  0.9],
           "LTR_ratio":[ 0.1 ,  0.15,  0.2 ,  0.25,  0.3 ,  0.35,  0.4 ,  0.45,  0.5 , 0.55,  0.6 ,  0.7,  0.8,  0.9]
           }

#[ 'upstream_intron_length', 'Exon_start_cons',  'Exon_end_cons', 'score_5p', 'score_3p']

use_log = ['upstream_intron_length', 'downstream_intron_length', 'Upstream_unaligned_ESS ', 'Upstream_aligned_ESS ', 'Upstream_aligned_RescueESS ', 'Downstream_unaligned_ESS ', 'Downstream_aligned_ESS ', 'Downstream_aligned_RescueESS ']


stats_numbers = {"Rs_score":[]}




for i in stats_fields:
    
    stats_numbers[i] = []



    
splice_site_scores = {}
if os.path.exists(maxent_file):
    
    print "Using maxent file! running second step!"
    maxent_here = True
    for line in open(maxent_file):
        
        
        ll = line.strip().split("\t")
    
        #printll
        
        
        splice_site_scores[ll[0]] = ll[1]
            #raw_input()
        



# input table is produced in R 
#

k562_skip = {}

control_analysis_pairs = {}

run_all =True

if run_all:
    
    print "running All!"
    with open(input_table) as csvfile:
        reader = csv.DictReader(csvfile)
        
        fields = reader.fieldnames
        
        '''print fields,"\n"
        print fields[:28],"\n"
        
        print fields[29:238],"\n"
        
        
        print fields.index("UPF2_ENCSR810FHY")
        
        raw_input()'''
        
        start = 0
        
        end = 0
        for i in fields:
    
            found_sample = False
    
            if (i.find("ENCS") > -1) and (i.find("_") > -1):
    
                found_sample = True
    
            if found_sample and start == 0 :
    
                start = fields.index(i)
    
            elif found_sample == False and start !=0:
    
                end = fields.index(i)
                break
    
        print "start",start
    
        print "end",end
    
        print fields[start:end]        
    
        #raw_input()
    
        for field in fields[start:end]:
            
            if field.find("_") == -1:
                        continue
            
            f_S = field.split("_")
            
            if len(f_S) == 2:
                print field
                sample_,experiment = f_S
            
            
            
            
            if sample_ == "Non": 
            
                print "sample is non"
                continue            
    
            if experiment not in control_id:
                print experiment
                
                print "Cant find control! ",field
                
                #raw_input("pause")
                
            else:
                
                samples.append(field)
                
                data_[field] = {}
                #print "field",field
                #print "experiment control is :",control_id[experiment]
                
                #raw_input()
                
                
        
        
        #sys.exit()           
    
        for row in reader:
            
            #print(row['first_name'], row['last_name'])
            count+=1
            
            id_ = row["skip_id"]
            
            ######## filtering low count data!!!
            if float(row["max_counts"]) < minimum_Reads:
                continue
            
            
            
            
            stats_dict = {}
            
            for stat_ in stats_fields:
                
                
                if stat_ != "Rs_score":
                    #print stat_
                    #print row[stat_]
                    #print row
                    
                    if row[stat_] == "NA":
                        continue
                    
                    elif row[stat_].strip() == "":
                        #stats_dict[stat_] = 0.0
                        continue
                    else:
                        
                        if stat_ in use_log:
                            stats_dict[stat_] = numpy.log(float(row[stat_]))
                                                
                            stats_numbers[stat_].append(numpy.log(float(row[stat_])))                    
                        
                        else:    
                            stats_dict[stat_] = float(row[stat_])
                            
                            stats_numbers[stat_].append(float(row[stat_]))
                else:     
            
            
                    rs_score = float(splice_site_scores[row["Rs_site_seq"]])
                    #rs_score = 
                    
                    stats_dict["Rs_score"] = rs_score
                    stats_numbers["Rs_score"].append(rs_score)
            
            
            for sam_ in samples:
                
                #print "sam_",sam_
                
                if len(sam_.split("_")) == 3:
                    
                    sample_,experiment = sam_.split("_")[1:]
                else:
                    sample_,experiment = sam_.split("_")
                
                
                #print control_id
                
                control_sample = "Non_"+control_id[experiment]
                
                
                '''
                print sample_,experiment
                
                #print control_analysis_pairs[experiment]
                
                print "control_sample",control_sample
                print float(row[control_sample])
                print float(row[sam_])
                
                #raw_input()
                
                
                #'''
                
                
                #data_[sam_][row["skip_id"]] = [float(row["Rs_site_score"]),float(row[control_sample]) - float(row[sam_]),row["psi_constitutive"]]
                
                control_psi = float(row[control_sample])
                
                delta_psi = float(row[control_sample]) - float(row[sam_])
                
                
                
                psi_cons = ""
                
                if control_psi > 0.98:
                    
                    psi_cons = "Cons"
                elif control_psi < 0.98 and control_psi > 0.9:
                    
                    psi_cons = "SemiCons"
                    
                else:
                    psi_cons = "Alt"
                    
                    
                
                
                #print  [row['eoi_chr'],row['eoi_start'],row['eoi_end'],row['x'],row['eoi_strand']]
            
                #raw_input()
                # 'eoi_chr', 'eoi_start', 'eoi_end', 'x', 'eoi_strand'
                #k562_skip[row["skip_id"]] =  [row['eoi_chr'],row['eoi_start'],row['eoi_end'],row['x'],rs_score,row['eoi_strand']]            
                    
                #difference_ratio = float(row["Rs_site_score"])/float(row["ALT_score"])
    
                data_[sam_][row["skip_id"]] = [stats_dict,delta_psi,psi_cons,row["max_counts"]]
                
                #full_data_csv.writerow([sam_,row["skip_id"],rs_score,delta_psi,psi_cons,row["sum_counts"]])
    
    
    # {"range":0,"limits":0,"step":0}
    
    
    msgpack.dump(stats_numbers,open("%s_stats_numbers_full.p"%(name),"w"))
    
    
    msgpack.dump(data_,open("%s_data_full.p"%(name),"w"))
    
    print "Dumped data dictionary"
    #pickle.dump(k562_skip,open("k562_skip.p","w"))
    
    '''import msgpack,csv,numpy
    
    
    data = msgpack.load(open("ENCODE_skip_V1_stats_numbers_full.p"))
    '''
else:
    stats_numbers =  msgpack.load(open("%s_stats_numbers_full.p"%(name)))



highlysig = {}

stats_details = {}
for stat_ in stats_numbers:
    
    print stat_
    
    if stat_ not in stats_fields:
        continue
    
    for i in stats_numbers[stat_]:
        
        if numpy.isinf(i) or i == numpy.Inf or i=="inf":
            
            #print "not finite!",i
            stats_numbers[stat_].remove(i)
    
    
    '''if set(stats_numbers[stat_]) < 20:
        stat_step.append(numpy.percentile(stats_numbers[stat_],percentile))
        
    '''
    stat_step = []    
    
    
    if stat_ in ranges_:
        
        stat_step = ranges_[stat_]
        
    else:
        for percentile in range(0,100,5):
            
            #numpy.percentile(k,95)
            stat_step.append(numpy.percentile(stats_numbers[stat_],percentile))
            
            if numpy.isnan(numpy.percentile(stats_numbers[stat_],percentile)):
                
                break
            # stat_step.append(numpy.percentile(stats_numbers[stat_],percentile))
            
        # stats_numbers.keys()
        # numpy.percentile(stats_numbers['Exon_end_cons'],5)
        #print set(stats_numbers[stat_])
        
        
        
        
        if len(set(stat_step)) < 8:
            
            stat_step = []
            
            if len(set(stats_numbers[stat_])) < 10:
                
                min_ = round(min(stats_numbers[stat_]))
                            
                max_ = math.ceil(max(stats_numbers[stat_]))
                
                step_ = 0.5         
                
                for i in numpy.arange(min_,max_,step_):
                                
                    print i 
                    
                    stat_step.append(i)
                
            else:
                
                
                
                min_ = round(min(stats_numbers[stat_]))
                
                max_ = math.ceil(max(stats_numbers[stat_]))
                
                step_ = float(max_ - min_)/10
                
                for i in numpy.arange(min_,max_,step_):
                    
                    print i 
                    
                    stat_step.append(i)
                    
            
        
    stats_details[stat_] = stat_step
    
    print stat_
    
    print stat_step
    
    print max(stats_numbers[stat_])
    print min(stats_numbers[stat_])
    print numpy.average(stats_numbers[stat_])
    #raw_input()




#sys.exit("Only generating full data pickle!")

            
### out files with data on rs scores
rs_Data_csv = csv.writer(open("%s_stats1_rsdata.csv"%(name),"w"))

highlysig_data = csv.writer(open("%s_highlysig_rsdata.csv"%(name),"w"))
     
        
### out files with data on psi data 
psi_Data_csv = csv.writer(open("%s_stats1_psidata.csv"%(name),"w"))

data_ =  msgpack.load(open("%s_data_full.p"%(name)))



write_raw_dotplots = True


'''if write_raw_dotplots:
    
    for sample_ in data_: 
        
        full_data_csv = csv.writer(open("%s_stats1_fulldata.csv"%(sample_),"w"))
        
        for skipid in data_[sample_]:
            
            full_data_csv.writerow(data_[sample_][skipid])
        [sam_][row["skip_id"]] = [stats_dict,delta_psi,psi_cons,row["sum_counts"]]
                



'''

rs_data = []

for stat_s in stats_numbers:

    if stat_s not in stats_fields:
        
        print "skipping ",stat_s
        continue
    
    
    print "Analysing ",stat_s 


    #continue
    for cons_var in ["Cons","Alt"]:
        
        
        
        psi_data = []
                    
        for Sample_ in data_:
            
            
            #print "Analysing: ",Sample_
                        
            #print "total data points: ",len(data_[Sample_])
            
            #print count
            
            ############## find best RS cut off
            
            best_rs_threshold = 0
            best_rs_threshold_tval = 0
            best_rs_threshold_pval = 2
            
            
            
            
            
            # for each recursive score threshold, make lists of data and calculate top vs bottom using t test, final all pvalues and statistics
            #print round(stats_details[stat_]["min"]),math.ceil(stats_details[stat_]["max"]),stats_details[stat_]["step"]
            for rs_thresh in  stats_details[stat_s]:
                #for rs_thresh in  range(-20,10,1):
                
                top_data = []
                bottom_data = []
                
                
                #print "rs_thresh",rs_thresh
                
                for val in data_[Sample_]:
                    
                    [stats_dict,psi_EIF4A3,psi_constitutive,max_count] = data_[Sample_][val]
                    
                    
                    try:
                        Rs_site_score = stats_dict[stat_s]
                    
                    except:
                    
                        continue
                    '''
                    print type(Rs_site_score)
                    print type(psi_EIF4A3)
                    
                    print Rs_site_score
                    
                    raw_input()'''
                    
                    if psi_constitutive != cons_var:
                        continue
                    
                    
                    '''if max_count < 15:
                        print "skipping dis"
                        continue
                    '''
                
                    if Rs_site_score > rs_thresh:
                        top_data.append(psi_EIF4A3)
                    else:
                        bottom_data.append(psi_EIF4A3)
                        
                        
                '''
                print rs_thresh
                print len(top_data)
                print len(bottom_data)
                
               
                #print top_data
                
                #print bottom_data
                
                #raw_input()#'''
                
                t, p = stats.ttest_ind(top_data, bottom_data, equal_var=False)
                
                
                wilk_s,wilk_p = stats.ranksums(top_data, bottom_data)
                
                
                rs_Data_csv.writerow([Sample_,stat_s,cons_var,rs_thresh,t,p,wilk_s,wilk_p,len(top_data),len(bottom_data)])
                
                
                if stat_s not in highlysig:
                    
                    highlysig[stat_s] = []
                    
                    
                    
                    
                # create a seperate file for those who are highly significant    
                if abs(t) >= 6.7 and (Sample_ not in highlysig[stat_s]):
                    
                    highlysig[stat_s].append(Sample_)
                    highlysig_data.writerow([Sample_,stat_s,cons_var,rs_thresh,t,p,len(top_data),len(bottom_data)])
                
  
    
       