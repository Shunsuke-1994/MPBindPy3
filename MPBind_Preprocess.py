import subprocess, re, os, sys, logging, time, datetime

################################### Read and Initialization Parameters #####################
# -Infile: Input file
# -t: infile type (Fastq or Fasta)
# -Forward_primer
# -Reverse_primer
# -primer_max_mismatch
# -Outfile

def function_print_usage():
    print('Usage:'+'\n')
    print('python MPBind_Preprocess.py -Infile -t <FASTA or FASTQ> -Forward_primer <Forward Primer sequence> -Reverse_primer <Reverse Primer sequence> -primer_max_mismatch <max mismatch allowed> -Outfile')
    
# Default Parameters
w_parameter = {}
w_parameter['-Infile'] = ''
w_parameter['-t'] = ''
w_parameter['-Forward_primer'] = ''
w_parameter['-Reverse_primer'] = ''
w_parameter['-primer_max_mismatch'] = ''
w_parameter['-Outfile'] = ''

Parameter_list = sys.argv[1:]

for i in range(len(Parameter_list)-1):
    if(Parameter_list[i] in w_parameter):
        w_parameter[Parameter_list[i]] = Parameter_list[i+1]
    else:
        pass

# Check Parameters
print(w_parameter)
if('' not in w_parameter.values()):
    pass
else:
    function_print_usage()
    sys.exit()

scriptPath = os.path.abspath(os.path.dirname(__file__))  ## absolute script path

print('Initialization Parameters Completed ...')

#########################################
if(w_parameter['-t'].upper() == 'FASTQ'):
    cmd = 'python'+'\t'+scriptPath+'/bin/From_Fastq_to_plain_text_and_transform_to_sense_seq.py'+'\t'+w_parameter['-Infile']+'\t'+w_parameter['-Forward_primer']+'\t'+w_parameter['-Reverse_primer']+'\t'+w_parameter['-primer_max_mismatch']+'\t'+w_parameter['-Outfile']
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if process.returncode != 0:
        print("Error executing command:")
        print(process.stderr)
        sys.exit(1)
elif(w_parameter['-t'].upper() == 'FASTA'):
    cmd = 'python'+'\t'+scriptPath+'/bin/From_Fasta_to_plain_text_and_transform_to_sense_seq.py'+'\t'+w_parameter['-Infile']+'\t'+w_parameter['-Forward_primer']+'\t'+w_parameter['-Reverse_primer']+'\t'+w_parameter['-primer_max_mismatch']+'\t'+w_parameter['-Outfile']
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if process.returncode != 0:
        print("Error executing command:")
        print(process.stderr)
        sys.exit(1)
else:
    print('Parameter Error ...')
    function_print_usage()
    sys.exit()