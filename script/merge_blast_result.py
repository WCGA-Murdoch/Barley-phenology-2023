# This script reads the BLAST result from files and write them to a merge file
# Require the Bio packaged to read the fasta file
# In the same directory of the script, there are 2 folder call
#    sequences: contains protein sequecens that used to blast against the database
#    blastoutput: blast result files


import os,re
from Bio import SeqIO

# Get the current directory of the script
base_dir = os.path.dirname(os.path.abspath(__file__))

#Specify the theshhold
len_thresh = 0.8 # match length need to be > 80%
iden_thresh = 0.95

# Append the path
sequence_dir= os.path.join(base_dir,'sequences')
bastoutput_dir = os.path.join(base_dir,'blastoutput')
output_file = os.path.join(base_dir,'blast_result_merged.txt')
# Write a blank file as a place holder
with open(output_file,'w') as o:  o.write('')


# Get the original protein sequence list and sort
sequences = os.listdir(sequence_dir)
sequences.sort()


### Merge the blast result
for sequence_file in sequences:
    # Get the gene id and gene length from the fasta sequence
    writeout = []
    sequence_path = os.path.join(sequence_dir,sequence_file)
    fasta_sequences = SeqIO.parse(open(sequence_path),'fasta')
    fasta = list(fasta_sequences)[0]
    gene_id, sequence_string = fasta.id, str(fasta.seq)
    gene_len = len(sequence_string)
    CHR = re.findall('.[1-7]H',gene_id)[0].replace('.','')
    writeout = [gene_id,str(gene_len)]

    # Get the blast result from the individual result file
    blast_result_file = os.path.join(bastoutput_dir,sequence_file+'.txt')
    hitfound = False
    # Check if there is hit found from blast
    # Bast fmt6 output 'Query','Reference','Identity','Length','Mismatches','Gaps','qstart','qend','sstart','send','E','bitscore'
    if os.path.exists(blast_result_file):
        if os.path.getsize(blast_result_file) > 0:
            with open(blast_result_file,'r') as f:
                for line in f:
                    if "#" not in line:
                        line = str(line)
                        line = re.sub('\n','',line)
                        line = re.split('\t',line)
                        # Try to numerize blast length and identity
                        try:
                            blast_len = int(line[3])
                            blast_iden  = float(line[2])
                        except:
                            blast_len = 0
                            blast_iden = 0
                        if blast_len >= gene_len*len_thresh and blast_iden >= iden_thresh*100:
                            blastgeneid = re.split('\.',line[1])[0]
                            identity = str(blast_iden)
                            coverage = str(blast_len)
                            if 'contig' not in blastgeneid: #If the result not in corrected conitig
                                # Then match the chromosome
                                CHR_BLAST = re.findall('chr[1-7]H',blastgeneid)[0].replace('chr','')
                                if CHR_BLAST == CHR:
                                    if blastgeneid not in writeout: # just incase of duplicates
                                        writeout += [blastgeneid,identity,coverage]
                            else:
                                if blastgeneid not in writeout:
                                    writeout += [blastgeneid,identity,coverage]

    # Write the result [Query GeneID,Gene length,Blast geneID, Blast identity, Blast length...]
    if len(writeout) > 2: # Onlywrite matched result
        with open(output_file,'a') as o:
            o.write('\t'.join(writeout)+'\n')


# Consolidate the blast table into defined columns
# [GeneID_GsRTD,GeneID_Variety,GeneID_Morex,Matchscore,LinearPanID,Source']
cleaned_outputfile = os.path.join(base_dir,'blast_result_cleaned.txt')

# Write empty files as place holder
with open(cleaned_outputfile,'w') as f:   f.write('')

gene_list = []
geneid_mapped = []

with open(output_file,'r') as f:
    for line in f:
        line = re.split('\t',line)
        geneid_var = line[0][0:len(line[0])-5] # Truncate the gemoma result text
        l = len(line)
        if l > 2: # At least 1 gene matched
            for i in range(int((l-2)/3)):
                geneid = line[i*3+2]
                if geneid not in geneid_mapped:
                    geneid_mapped +=[geneid]
                    identity  = str(line[i*3+3])
                    with open(cleaned_outputfile,'a') as g:
                        geneid_morex = re.split('_',geneid_var)[1] # Truncate the variety name
                        writeout = [geneid, geneid_var,identity]
                        g.write('\t'.join(writeout)+'\n')
