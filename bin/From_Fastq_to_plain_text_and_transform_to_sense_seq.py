import subprocess, sys    # subprocess.getstatusoutput(cmd)

def function_reverse_complement(seq):
    seq=seq.upper()
    swap_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    newword=''.join(swap_dict[letter] for letter in reversed(seq))
    return newword

#
Fastq_file=sys.argv[1]
Forward_primer=sys.argv[2].upper()
Reverse_primer=sys.argv[3].upper()
primer_max_mismatch=int(sys.argv[4])
ouput_file=sys.argv[5]

#
infile=open(Fastq_file,'r')
line_name=infile.readline()
line_seq=infile.readline().rstrip()
Random_Region_len=len(line_seq)-len(Forward_primer)-len(Reverse_primer)
infile.close()

####################################
# len(seq_1)==len(seq_2)
# return Yes or No (mismatch<=cutoff)
def function_less_or_eq_mismatch_logical(seq_1,seq_2,mismatch_cutoff):
    if(len(seq_1)!=len(seq_2)):
        print('Error... Len_1!=Len2')
    else:
        pass
    logical_less_or_eq_mismatch='Yes'
    N_mismatch=int(0)
    seq_len=len(seq_1)
    for i in range(seq_len):
        if(N_mismatch<=mismatch_cutoff):
            if(seq_1[i]!=seq_2[i]):
                N_mismatch=N_mismatch+1
            else:
                pass
        else:
            logical_less_or_eq_mismatch='No'
            break
    if(N_mismatch>mismatch_cutoff):
        logical_less_or_eq_mismatch='No'
    else:
        pass
    return logical_less_or_eq_mismatch

infile=open(Fastq_file,'r')
outfile=open(ouput_file,'w')
while True:
    line_name=infile.readline()
    if not line_name:
        break
    line_seq=infile.readline()
    line_temp=infile.readline()
    line_Q=infile.readline()
    
    seq_whole=line_seq.rstrip().upper()
    if(seq_whole.find('N')==(-1) and seq_whole.find('.')==(-1)):
        primer_region=seq_whole[0:len(Forward_primer)]
        random_region=seq_whole[len(Forward_primer):(len(Forward_primer)+Random_Region_len)]
        if(function_less_or_eq_mismatch_logical(primer_region,Forward_primer,primer_max_mismatch)=='Yes'):
            outfile.write(random_region.upper()+'\n')
        else:
            primer_region=seq_whole[0:len(Reverse_primer)]
            random_region=seq_whole[len(Reverse_primer):(len(Reverse_primer)+Random_Region_len)]
            if(function_less_or_eq_mismatch_logical(primer_region,Reverse_primer,primer_max_mismatch)=='Yes'):
                outfile.write(function_reverse_complement(random_region.upper())+'\n')
            else:
                pass
    else:
        pass

infile.close()
outfile.close()


