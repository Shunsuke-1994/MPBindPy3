import subprocess, sys, os

scriptPath=os.path.abspath(os.path.dirname(__file__));  ## absolute script path

All_files=sys.argv[1].split(',')
outfile_name=sys.argv[2]

outfile=open(outfile_name,'w')
header='Motif'+'\t'+'Fraction_list'+'\t'+'P_value[1-sided]'+'\t'+'Z_Score'+'\n'
outfile.write(header)

round_list=list(range(len(All_files)))

w_unique_motif={}
w_complex={}

for i in round_list:
    infile=open(All_files[i],'r')
    header_temp=infile.readline()
    for line in infile:
        motif=line.split()[0]
        Fraction=line.split()[3] # *** by_substring ***
        w_unique_motif[motif]=''
        w_complex[motif+':'+str(i)]=Fraction
    infile.close()

for k in w_unique_motif.keys():
    outfile.write(k+'\t')
    fraction_list=[]
    for i in round_list:
        fraction_list.append(w_complex[k+':'+str(i)])
    outfile.write(','.join(fraction_list)+'\t')
    fraction_list=[str(j) for j in fraction_list]
    round_this_list=[str(j) for j in round_list]
    cmd='Rscript'+'\t'+scriptPath+'/Spearman_1_sided_P_value.R'+'\t'+','.join(fraction_list)+'\t'+','.join(round_this_list)
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    P_value=proc.stdout.strip()
    cmd='Rscript'+'\t'+scriptPath+'/from_P_value_to_Z_Score_command.R'+'\t'+str(P_value)
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    Z_Score=proc.stdout.strip()
    outfile.write(str(P_value)+'\t'+str(Z_Score)+'\n')

outfile.close()

