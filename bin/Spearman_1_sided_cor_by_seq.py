import commands,sys,os

scriptPath=os.path.abspath(os.path.dirname(__file__));  ## absolute script path

All_files=sys.argv[1].split(',')
outfile_name=sys.argv[2]

outfile=open(outfile_name,'w')
header='Motif'+'\t'+'Fraction_list'+'\t'+'P_value[1-sided]'+'\t'+'Z_Score'+'\n'
outfile.write(header)

round_list=range(len(All_files))

w_unique_motif={}
w_complex={}

for i in round_list:
    infile=open(All_files[i],'r')
    header_temp=infile.readline()
    line=infile.readline()
    while(line):
        motif=line.split()[0]
        Fraction=line.split()[6] # *** by_seq ***
        w_unique_motif[motif]=''
        w_complex[motif+':'+str(i)]=Fraction
        line=infile.readline()
    infile.close()

for k in w_unique_motif.keys():
    outfile.write(k+'\t')
    fraction_list=[]
    for i in round_list:
        fraction_list.append(w_complex[k+':'+str(i)])
    outfile.write(','.join(fraction_list)+'\t')
    fraction_list=[str(j) for j in fraction_list]
    round_this_list=[str(j) for j in round_list]
    cmd='python'+'\t'+scriptPath+'/Spearman_1_sided_P_value.py'+'\t'+','.join(fraction_list)+'\t'+','.join(round_this_list)
    a=commands.getstatusoutput(cmd)
    P_value=a[1]
    cmd='python'+'\t'+scriptPath+'/from_P_value_to_Z_Score_command.py'+'\t'+str(P_value)
    a=commands.getstatusoutput(cmd)
    Z_Score=a[1]
    outfile.write(str(P_value)+'\t'+str(Z_Score)+'\n')

outfile.close()

