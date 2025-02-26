# Usage: python nmer_frequency_single_version.py <nmer> <infile> <outfile>
# E.g, python nmer_frequency_single_version.py 6 Round_0_unambiguous_SELEX_Seq_uppercase.txt 1_frequency_6mer_Round_0_SELEX.txt

import sys

if(len(sys.argv)!=4):
    print('Usage: python nmer_frequency_single_version.py <nmer> <infile> <outfile>')
    sys.exit(1)
else:
    pass

user_defined_nmer=int(sys.argv[1])
infile_SELEX_Seq=sys.argv[2]
outfile_freq=sys.argv[3]

######## Function for Generate nmers #############################
def product(*args, **kwds):
    # product('ABCD', 'xy') --> Ax Ay Bx By Cx Cy Dx Dy
    # product(range(2), repeat=3) --> 000 001 010 011 100 101 110 111
    pools = list(map(tuple, args)) * kwds.get('repeat', 1)
    result = [[]]
    for pool in pools:
        result = [x+[y] for x in result for y in pool]
    for prod in result:
        yield tuple(prod)

keys_tuple_nmer=list(product('ATCG', repeat=user_defined_nmer))

def function_parse_nmer(keys_tuple_nmer):
    nmer_list_this=[]
    for i in keys_tuple_nmer:
        nmer_list_this.append(''.join(i))
    return nmer_list_this

nmer_list=function_parse_nmer(keys_tuple_nmer)  # ***
motif_len=int(user_defined_nmer)

def Function_check_ATCG(read,w_ATCG):
    for i in range(len(read)):
        if(read[i] in w_ATCG):
            pass
        else:
            print('Error: Non-{ATCG} charactor found in the sequence')
            print(read)
            exit(1)

w_ATCG={}
w_ATCG['A']=''
w_ATCG['T']=''
w_ATCG['C']=''
w_ATCG['G']=''

# print('motif_len: '+str(motif_len))

#############################################################
w_nmer_by_substring={}
for i in nmer_list:
    w_nmer_by_substring[i]=int(0)

w_nmer_by_seq={}
for i in nmer_list:
    w_nmer_by_seq[i]=int(0)

##########################################
def function_split_into_substring(substring_len,seq):
    substring_list=[]
    for i in range(len(seq)-int(substring_len)+1):
        substring_list.append(seq[i:(i+int(substring_len))])
    return substring_list
#########################################
sum_all_seq_counts=int(0)

infile=open(infile_SELEX_Seq,'r')
# header_temp=infile.readline()

for line in infile:
    seq=line.split()[0].upper()
    Function_check_ATCG(seq,w_ATCG)
    substring_list_this_seq=function_split_into_substring(motif_len,seq)
    for k in substring_list_this_seq:
        w_nmer_by_substring[k]=w_nmer_by_substring[k]+1

    w_unique={}
    for kk in substring_list_this_seq:
        w_unique[kk]=''
    for m in w_unique.keys():
        w_nmer_by_seq[m]=w_nmer_by_seq[m]+1

    sum_all_seq_counts=sum_all_seq_counts+1

infile.close()

#### Outfile ###############
sum_all_substring_counts=int(0)

for i in nmer_list:
    sum_all_substring_counts=sum_all_substring_counts+w_nmer_by_substring[i]

outfile=open(outfile_freq,'w')
header='Motif'+'\t'+'Counts_by_substring'+'\t'+'Total_counts_all_substring'+'\t'+'Fraction_by_substring'+'\t'+'Counts_by_reads'+'\t'+'Total_reads'+'\t'+'Fraction_by_reads'+'\n'

outfile.write(header)

for k in nmer_list:
    outfile.write(k+'\t')
    Counts_by_substring=w_nmer_by_substring[k]
    Fraction_by_substring=float(Counts_by_substring)/float(sum_all_substring_counts)
    outfile.write(str(Counts_by_substring)+'\t'+str(sum_all_substring_counts)+'\t'+str(Fraction_by_substring)+'\t')
    Counts_by_reads=w_nmer_by_seq[k]
    Fraction_by_reads=float(Counts_by_reads)/float(sum_all_seq_counts)
    outfile.write(str(Counts_by_reads)+'\t'+str(sum_all_seq_counts)+'\t'+str(Fraction_by_reads)+'\n')

outfile.close()













