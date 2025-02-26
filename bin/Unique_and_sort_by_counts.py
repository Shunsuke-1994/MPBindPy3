import subprocess, sys    # subprocess.run() for commands

def function_unique_and_sort_by_counts(infile_name,outfile_name):
    w_unique={}
    infile=open(infile_name,'r')
    for line in infile:
        seq=line.rstrip()
        if(seq in w_unique):
            w_unique[seq]=w_unique[seq]+1
        else:
            w_unique[seq]=int(1)
    infile.close()
    all_list=[]
    for i in w_unique.keys():
        all_list.append([w_unique[i],i])
    all_list.sort()
    all_list.reverse()
    outfile=open(outfile_name,'w')
#   header='Unique_Sense_Seq'+'\t'+'Counts'+'\n'
#   outfile.write(header)
    for i,j in all_list:
        outfile.write(j+'\t'+str(i)+'\n')
    outfile.close()

##############
infile_name=sys.argv[1]
outfile_name=sys.argv[2]

function_unique_and_sort_by_counts(infile_name,outfile_name)







