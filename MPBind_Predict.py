import subprocess, re, os, sys, logging, time, datetime

start_time = time.time()
################################### Read and Initialization Parameters #####################
# -Train
# -Aptamer
# -Sort
# -Out (outfile folder)

def function_print_usage():
    print('Usage:'+'\n')
    print('python MPBind_Predict.py -Train <train file: output from MPBind_Predict> -Aptamer <input aptamer file> -Sort <TRUE or FALSE> -Out <outfile>'+'\n')

# Default Parameters
w_parameter = {}
w_parameter['-Train'] = ''
w_parameter['-Aptamer'] = ''
w_parameter['-Sort'] = 'FALSE'
w_parameter['-Out'] = ''

Parameter_list = sys.argv[1:]

#
for i in range(len(Parameter_list)-1):
    if(Parameter_list[i] in w_parameter):
        w_parameter[Parameter_list[i]] = Parameter_list[i+1]
    else:
        pass

# check Parameters
print('The parameters of MPBind (Prediction):')
for k in w_parameter:
    print(k+'='+w_parameter[k])
print('\n')

if('' not in w_parameter.values()):
    pass
else:
    print('Missing parameters:')
    for i in w_parameter.keys():
        if(w_parameter[i]==''):
            print(i)
    print('*****************************')
    function_print_usage()
    print('*****************************')
    sys.exit()

if(w_parameter['-Sort'].upper()=='TRUE' or w_parameter['-Sort'].upper()=='FALSE'):
    pass
else:
    print('Unrecognized user input parameter for -Sort option')
    print('-Sort: TRUE or FALSE; Default=FALSE')
    sys.exit()

scriptPath = os.path.abspath(os.path.dirname(__file__))  ## absolute script path

print('Initialization Parameters Complete ...')

def function_Run_time(start_time):
    runningTime = time.time() - start_time
    print("Time Elapsed:  %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60))

def Function_Checking_Running_Status(returncode, output):
    if(returncode == 0):
        print('Running Status Checking: PASS ...')
    else:
        print('Encountering a problem:')
        print(output)
        exit(1)

####################################################################################

# Determine Z1_Z2_Z3_Z4 or Z3_Z4

with open(w_parameter['-Train'], 'r') as infile:
    header = infile.readline()

# print header.split()

if(len(header.split()) == 6):
    print('Meta-Z-Score: Combined 4 kinds of Z-Scores (Z1, Z2, Z3 and Z4)')
    cmd = 'python'+'\t'+scriptPath+'/bin/Aptamer_Scan_Combined_Z1_Z2_Z3_Z4.py'+'\t'+w_parameter['-Train']+'\t'+w_parameter['-Aptamer']+'\t'+w_parameter['-Sort']+'\t'+w_parameter['-Out']
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    function_Run_time(start_time)
    Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)    
else:
    print('The user did not provide final Control-Seq for training. Thus MPBind calculated the Meta-Z-Score using (combined Z3 and Z4) and skipped Z1 and Z2')
    print('Meta-Z-Score: Combined 2 kinds of Z-Scores (Z3 and Z4)')
    cmd = 'python'+'\t'+scriptPath+'/bin/Aptamer_Scan_Combined_Z3_Z4.py'+'\t'+w_parameter['-Train']+'\t'+w_parameter['-Aptamer']+'\t'+w_parameter['-Sort']+'\t'+w_parameter['-Out']
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    function_Run_time(start_time)
    Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)

print('Prediction Process is finished ...')
runningTime = time.time() - start_time
print("Total Running Time:  %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60))