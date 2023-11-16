# This script reads cleaned merged BLAST results from a table and match them with the TPM by geneid
# The script output a file name tpm_output
# The TPM files are put in a folder called tpm and in the same place as the script

import os,re
base_dir = os.path.dirname(os.path.abspath(__file__))

cleaned_outputfile = os.path.join(base_dir,'blast_result_cleaned.txt')
tpmdir = os.path.join(base_dir,tpm)
output_tpm = os.path.join(base_dir,'tpm_output.txt')

# Create a dict with tpm name for matching
tpmval_template = dict.fromkeys(
            [
            'Geneid_GsRTDs',
            'MorexV3_Variety',
            'MorexV3_Geneid',
            'Ca1',
            'Ca2',
            'Ca3',
            'Co1',
            'Co2',
            'Co3',
            'In1',
            'In2',
            'In3',
            'Ro1',
            'Ro2',
            'Ro3',
            'Sh1',
            'Sh2',
            'Sh3'
            ]
            ,''
        )
# Write out file with column names
for idx,value in tpmval_template.items():
    writeout.append(idx)
with open(outputf,'w') as f:
    f.write('\t'.join(writeout) + '\n')

# Iterate the cleaned blast result by line
with open(cleaned_outputfile,'r') as file:
    for line in file:
        # Create a dict from template
        tpmval = tpmval_template
        writeout = []
        if 'GeneID' not in line: # Exclude the first line
            line = re.split('\t',line)
            geneid = line[0]
            morex_varid = line[1]
            morexid = line[2]
            # Get the current variety name from the geneid and read the df
            current_var = re.split('_',geneid)[0]
            tpmfile = os.path.join(tpmdir,'Genes_TPM_' + current_var +'.csv')
            df = pd.read_csv(tpmfile)

            # Iterate through the tpm template
            for idx,value in tpmval.items():
                # Column name is varity.tissue
                colname = current_var+"."+idx
                # Get the gene ids
                if idx == 'Geneid_GsRTDs': tpmval[idx] = geneid
                elif idx == 'MorexV3_Variety': tpmval[idx] = morex_varid
                elif idx == 'MorexV3_Geneid': tpmval[idx] = morexid
                # Match the tpm
                else:
                    if colname in df.columns.tolist():
                        tpmval[idx] = str(df.loc[df["Unnamed: 0"]==geneid,colname].item())
                    else:
                        tpmval[idx] = 'NA'
            # Get the value into a list
            for idx,value in tpmval.items():
                writeout.append(tpmval[idx])
            # Write the output
            with open(outputf,'a') as f:
                f.write('\t'.join(writeout) + '\n')
