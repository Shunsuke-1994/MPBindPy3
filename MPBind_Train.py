import subprocess, re, os, sys, logging, time, datetime

start_time = time.time()

print('MPBind Training Starts ...')
print('\n')

Project_prefix = 'MPBind'
################################### Read and Initialization Parameters #####################
# -R0
# -RS  (SELEX rounds)
# -RC (control rounds)  -> 'NULL' or 'control_1,control_2,...'
# -nmer (5,6,7...)
# -U (user: '1': Unique reads only; '2': Redundant reads only; '3' both)
# -Out (outfile folder)

def function_print_usage():
    print('\n')
    print('Usage:'+'\n')
    print('python MPBind_Train.py -R0 <Round 0 library files> -RS <SELEX round files> -RC <Control round files> -nmer <motif_length> -U <1: unique reads only; 2: Redundant reads only; 3: Both> -Out <Outfile folder>'+'\n')
    print('Details:'+'\n')
    print('-R0 <Round 0 library files>'+'\t'+'e.g., Round_0.txt (plain text, each row contain one sequence, header=FALSE)')
    print('-RS  (SELEX rounds)'+'\t'+'e.g., Round_1.txt,Round_2.txt,Round_3.txt,... (plain text, each row contain one sequence, header=FALSE)')
    print('-RC <finial control round>'+'\t'+'e.g.,Control_3.txt (plain text, each row contain one sequence, header=FALSE)'+'--this is optional and default is NULL, if NULL, only Z3 and Z4 are calcuated and Meta-Z-Score will be combined by Z3 and Z4, skipping Z1 and Z2)')
    print('-nmer <motif_length>'+'\t'+'e.g, 5,6,7 <default is 6 (6mer)>')
    print('-U <1: Unique reads only; 2: Redundant reads only; 3 both> (default is 1)')
    print('-Out <Outfile folder, default is MPBind_Out>'+'\n')

# Parameter_list=sys.argv[1:]
def Function_Parse_User_Input(Parameter_list, w_default_parameter):
    Key_list = list(w_default_parameter.keys())  # ['-R0','-RS',...]
    w_user_parameter = {}
    for i in range(len(Parameter_list)):
        if(Parameter_list[i] in w_default_parameter):
            This_user_parameter = []
            for j in range(i+1, len(Parameter_list)):
                if((Parameter_list[j] in w_default_parameter)):   
                    w_user_parameter[Parameter_list[i]] = ''.join(This_user_parameter)
                    break
                elif(j == len(Parameter_list)-1):
                    This_user_parameter.append(Parameter_list[j])
                    w_user_parameter[Parameter_list[i]] = ''.join(This_user_parameter)
                    break
                else:
                    This_user_parameter.append(Parameter_list[j])
        else:
            pass

    for k in Key_list:
        if(k in w_user_parameter):
            pass
        else:
            w_user_parameter[k] = w_default_parameter[k]

    if('' not in w_user_parameter.values()):
        pass
    else:
        print('Please check the missing values [***]')
        for k in w_default_parameter.keys():
            if(w_user_parameter[k] != ''):
                print(k+'='+w_user_parameter[k])
            else:
                print(k+'='+'***')
        print('Incomplete parameters input or parameter input contains unrecognized characters!')
        print('For example, characters directly copied from PDF might cause this issue ...')
        print('Please check your input ...')
        function_print_usage()
        sys.exit()
    if(w_user_parameter['-Out'][-1] == '/'):
        w_user_parameter['-Out'] = ''.join(w_user_parameter['-Out'][:-1])
    else:
        pass
    print('The parameters of MPBind (Training):')
    for k in w_default_parameter.keys():
        print(k+'='+w_user_parameter[k])
    print('\n')
    return w_user_parameter

#
def function_Run_time(start_time):
    runningTime = time.time() - start_time
    print("Time Elapsed:  %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60))

#########################

# Default Parameters
w_default_parameter = {}
w_default_parameter['-R0'] = ''
w_default_parameter['-RS'] = ''
w_default_parameter['-RC'] = 'NULL'
w_default_parameter['-nmer'] = '6' 
w_default_parameter['-U'] = '1'
w_default_parameter['-Out'] = 'MPBind_Out'


Parameter_list = sys.argv[1:]

w_parameter = Function_Parse_User_Input(Parameter_list, w_default_parameter)

scriptPath = os.path.abspath(os.path.dirname(__file__))  ## absolute script path

print('Initialization Parameters Completed ...')
print('\n')

##
def Function_Checking_Running_Status(returncode, output):
    if(returncode == 0):
        print('Running Status Checking: PASS ...')
    else:
        print('Encountering a problem:')
        print(output)
        exit(1)

########################## Generate Unique Reads ####################
print('Reads preparing  ...')

cmd = 'mkdir -p ' + w_parameter['-Out']
process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
returncode = process.returncode
output = process.stdout + process.stderr

Function_Checking_Running_Status(returncode, output)

if(returncode != int(0)):
    print('Making folder Error ...')
    print('Folder Already Exist ...')
    exit(1)
else:
    pass

redundant_file_list = [w_parameter['-R0']]
unique_file_list = []

for i in w_parameter['-RS'].split(','):
    redundant_file_list.append(i)
if(w_parameter['-RC'] != 'NULL'):
    redundant_file_list.append(w_parameter['-RC'])
else:
    pass

for i in redundant_file_list:
    unique_file_list.append(w_parameter['-Out']+'/Unique_Reads/'+i.split('/')[-1].split('.txt')[0]+'.unique')   # i.split('.txt')[0] -> i.split('/')[-1].split('.txt')[0]

if(w_parameter['-U'] == '1'):
    cmd = 'mkdir -p ' + w_parameter['-Out']+'/Unique_Reads'
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    function_Run_time(start_time)
    Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
    for i in range(len(redundant_file_list)):
        this_redundant_file = redundant_file_list[i]
        this_unique_file = unique_file_list[i]
        cmd = 'python'+'\t'+scriptPath+'/bin/Unique_and_sort_by_counts.py'+'\t'+this_redundant_file+'\t'+this_unique_file
        process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        function_Run_time(start_time)
        Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
        if(process.returncode != int(0)):
            print('Error in Generating Unique Reads from Redundant Reads ...')
            print('The Error file (redundant reads):  ')
            print(this_redundant_file)
            print(process.stdout + process.stderr)
            exit(1)
        else:
            pass
elif(w_parameter['-U'] == '2'):
    cmd = 'mkdir -p ' + w_parameter['-Out']+'/Redundant_Reads'
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    function_Run_time(start_time)
    Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
elif(w_parameter['-U'] == '3'):
    cmd = 'mkdir'+'\t'+w_parameter['-Out']+'/Redundant_Reads'
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    function_Run_time(start_time)
    Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
    cmd = 'mkdir'+'\t'+w_parameter['-Out']+'/Unique_Reads'
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    function_Run_time(start_time)
    Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
    for i in range(len(redundant_file_list)):
        this_redundant_file = redundant_file_list[i]
        this_unique_file = unique_file_list[i]
        cmd = 'python'+'\t'+scriptPath+'/bin/Unique_and_sort_by_counts.py'+'\t'+this_redundant_file+'\t'+this_unique_file
        process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        function_Run_time(start_time)
        Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
        if(process.returncode != int(0)):
            print('Error in Generating Unique Reads from Redundant Reads ...')
            print('The Error file (redundant reads):  ')
            print(this_redundant_file)
            print(process.stdout + process.stderr)
            exit(1)
        else:
            pass
else:
    print('Error: -U parameter')

print('Reads preparing are completed ...')
function_Run_time(start_time)

##################### nmer ######################

#print redundant_file_list # ['Random.R0.sense_seq.txt', 'Random_R1_SELEX.sense_seq.txt', 'Random_R2_SELEX.sense_seq.txt', 'Random_R3_SELEX.sense_seq.txt', 'Random_R3_Control.sense_seq.txt']
#print unique_file_list # ['./Peng_PRAS_Out/Unique_Reads/Random.R0.sense_seq.txt.unique', './Peng_PRAS_Out/Unique_Reads/Random_R1_SELEX.sense_seq.txt.unique', './Peng_PRAS_Out/Unique_Reads/Random_R2_SELEX.sense_seq.txt.unique', './Peng_PRAS_Out/Unique_Reads/Random_R3_SELEX.sense_seq.txt.unique', './Peng_PRAS_Out/Unique_Reads/Random_R3_Control.sense_seq.txt.unique']

print('Starting n-mer calculating ...')

nmer_list = w_parameter['-nmer'].split(',')

if(w_parameter['-U'] == '1'):
    for i in nmer_list:
        for j in unique_file_list:
            infile = j
            outfile = '/'.join(j.split('/')[:-1])+'/'+j.split('/')[-1]+'.'+i+'mer' # Fix
            cmd = 'python'+'\t'+scriptPath+'/bin/Header_free_nmer_frequency_single_version.py'+'\t'+i+'\t'+infile+'\t'+outfile
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
elif(w_parameter['-U'] == '2'):
    for i in nmer_list:
        for j in redundant_file_list:
            infile = j
            outfile = w_parameter['-Out']+'/Redundant_Reads/'+j.split('/')[-1].split('.txt')[0]+'.'+i+'mer'       # Fix  .split('.txt')[0] -> .split('/')[-1].split('.txt')[0] 
            cmd = 'python'+'\t'+scriptPath+'/bin/Header_free_nmer_frequency_single_version.py'+'\t'+i+'\t'+infile+'\t'+outfile
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
elif(w_parameter['-U'] == '3'):
    for i in nmer_list:
        for j in unique_file_list:
            infile = j
            outfile = '/'.join(j.split('/')[:-1])+'/'+j.split('/')[-1]+'.'+i+'mer'   # Fix
            cmd = 'python'+'\t'+scriptPath+'/bin/Header_free_nmer_frequency_single_version.py'+'\t'+i+'\t'+infile+'\t'+outfile
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
    for i in nmer_list:
        for j in redundant_file_list:
            infile = j
            outfile = w_parameter['-Out']+'/Redundant_Reads/'+j.split('/')[-1].split('.txt')[0]+'.'+i+'mer'   # Fix
            cmd = 'python'+'\t'+scriptPath+'/bin/Header_free_nmer_frequency_single_version.py'+'\t'+i+'\t'+infile+'\t'+outfile
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)

print('n-mer completed ...')

################ p-values and Z1,Z2,Z3,Z4  ###################

print('Starting p-values and Z-Scores calculating ...')

if(w_parameter['-U'] == '1'):
    folder_path = w_parameter['-Out']+'/Unique_Reads/'
    if(w_parameter['-RC'] != 'NULL'):
        for i in nmer_list:
            finial_SELEX = w_parameter['-RS'].split(',')[-1].split('/')[-1].split('.txt')[0]+'.unique.'+i+'mer'    # .split('.txt')[0] -> .split('/')[-1].split('.txt')[0] [Fix]
            finial_Control = w_parameter['-RC'].split('/')[-1].split('.txt')[0]+'.unique.'+i+'mer' # Fix
            Round_list = []
            for k in w_parameter['-RS'].split(','):
                Round_list.append(folder_path+k.split('/')[-1].split('.txt')[0]+'.unique.'+i+'mer') # Fix
            cmd = 'python'+'\t'+scriptPath+'/bin/Fisher_test_1_sided_SELEX_Vs_Control.py'+'\t'+folder_path+finial_SELEX+'\t'+folder_path+finial_Control+'\t'+folder_path+'Z1_Z2.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_substring.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z3.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_seq.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z4.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
    else:
        for i in nmer_list:
            Round_list = []
            for k in w_parameter['-RS'].split(','):
                Round_list.append(folder_path+k.split('/')[-1].split('.txt')[0]+'.unique.'+i+'mer') # Fix
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_substring.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z3.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_seq.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z4.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)

elif(w_parameter['-U'] == '2'):
    folder_path = w_parameter['-Out']+'/Redundant_Reads/'
    if(w_parameter['-RC'] != 'NULL'):
        for i in nmer_list:
            finial_SELEX = w_parameter['-RS'].split(',')[-1].split('/')[-1].split('.txt')[0]+'.'+i+'mer'  # Fix
            finial_Control = w_parameter['-RC'].split('/')[-1].split('.txt')[0]+'.'+i+'mer'  # Fix
            Round_list = []
            for k in w_parameter['-RS'].split(','):
                Round_list.append(folder_path+k.split('/')[-1].split('.txt')[0]+'.'+i+'mer') # Fix
            cmd = 'python'+'\t'+scriptPath+'/bin/Fisher_test_1_sided_SELEX_Vs_Control.py'+'\t'+folder_path+finial_SELEX+'\t'+folder_path+finial_Control+'\t'+folder_path+'Z1_Z2.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_substring.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z3.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_seq.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z4.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
    else:
        for i in nmer_list:
            Round_list = []
            for k in w_parameter['-RS'].split(','):
                Round_list.append(folder_path+k.split('/')[-1].split('.txt')[0]+'.'+i+'mer')   # Fix
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_substring.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z3.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_seq.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z4.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)

elif(w_parameter['-U'] == '3'):
    folder_path = w_parameter['-Out']+'/Unique_Reads/'
    if(w_parameter['-RC'] != 'NULL'):
        for i in nmer_list:
            finial_SELEX = w_parameter['-RS'].split(',')[-1].split('/')[-1].split('.txt')[0]+'.unique.'+i+'mer'     # Fix
            finial_Control = w_parameter['-RC'].split('/')[-1].split('.txt')[0]+'.unique.'+i+'mer'     # Fix
            Round_list = []
            for k in w_parameter['-RS'].split(','):
                Round_list.append(folder_path+k.split('/')[-1].split('.txt')[0]+'.unique.'+i+'mer')       # Fix
            cmd = 'python'+'\t'+scriptPath+'/bin/Fisher_test_1_sided_SELEX_Vs_Control.py'+'\t'+folder_path+finial_SELEX+'\t'+folder_path+finial_Control+'\t'+folder_path+'Z1_Z2.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_substring.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z3.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_seq.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z4.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
    else:
        for i in nmer_list:
            Round_list = []
            for k in w_parameter['-RS'].split(','):
                Round_list.append(folder_path+k.split('/')[-1].split('.txt')[0]+'.unique.'+i+'mer')   # Fix
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_substring.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z3.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_seq.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z4.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
    
    folder_path = w_parameter['-Out']+'/Redundant_Reads/'
    if(w_parameter['-RC'] != 'NULL'):
        for i in nmer_list:
            finial_SELEX = w_parameter['-RS'].split(',')[-1].split('/')[-1].split('.txt')[0]+'.'+i+'mer'        # Fix
            finial_Control = w_parameter['-RC'].split('/')[-1].split('.txt')[0]+'.'+i+'mer'                    # Fix
            Round_list = []
            for k in w_parameter['-RS'].split(','):
                Round_list.append(folder_path+k.split('/')[-1].split('.txt')[0]+'.'+i+'mer') # Fix
            cmd = 'python'+'\t'+scriptPath+'/bin/Fisher_test_1_sided_SELEX_Vs_Control.py'+'\t'+folder_path+finial_SELEX+'\t'+folder_path+finial_Control+'\t'+folder_path+'Z1_Z2.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_substring.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z3.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_seq.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z4.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
    else:
        for i in nmer_list:
            Round_list = []
            for k in w_parameter['-RS'].split(','):
                Round_list.append(folder_path+k.split('/')[-1].split('.txt')[0]+'.'+i+'mer')   # Fix
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_substring.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z3.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
            cmd = 'python'+'\t'+scriptPath+'/bin/Spearman_1_sided_cor_by_seq.py'+'\t'+','.join(Round_list)+'\t'+folder_path+'Z4.'+i+'mer'
            process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            function_Run_time(start_time)
            Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)

else:
    print('Error ... (-U)')

print('Individual Z-Scores for each n-mer are completed  ...')

############################ Combine Z-Score ##########################
if(w_parameter['-U'] == '1'):
    folder_path = w_parameter['-Out']+'/Unique_Reads/'
    cmd = 'python'+'\t'+scriptPath+'/bin/Folder_nmer_combined_Z_Score.py'+'\t'+folder_path+'\t'+w_parameter['-RC']+'\t'+w_parameter['-nmer']
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    function_Run_time(start_time)
    Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
elif(w_parameter['-U'] == '2'):
    folder_path = w_parameter['-Out']+'/Redundant_Reads/'
    cmd = 'python'+'\t'+scriptPath+'/bin/Folder_nmer_combined_Z_Score.py'+'\t'+folder_path+'\t'+w_parameter['-RC']+'\t'+w_parameter['-nmer']
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    function_Run_time(start_time)
    Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
elif(w_parameter['-U'] == '3'):
    folder_path = w_parameter['-Out']+'/Unique_Reads/'
    cmd = 'python'+'\t'+scriptPath+'/bin/Folder_nmer_combined_Z_Score.py'+'\t'+folder_path+'\t'+w_parameter['-RC']+'\t'+w_parameter['-nmer']
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    function_Run_time(start_time)
    Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
    folder_path = w_parameter['-Out']+'/Redundant_Reads/'
    cmd = 'python'+'\t'+scriptPath+'/bin/Folder_nmer_combined_Z_Score.py'+'\t'+folder_path+'\t'+w_parameter['-RC']+'\t'+w_parameter['-nmer']
    process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    function_Run_time(start_time)
    Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
else:
    print('Error ...combine Z-score')

print('Combined Z-score completed ...')

print('Compiling Training Result files ...')
################################### Copy ################################
outfile_list = []

if(w_parameter['-U'] == '1'):
    folder_path = w_parameter['-Out']+'/Unique_Reads/'
    for i in nmer_list:
        cmd = 'cp'+'\t'+folder_path+'Combined_Z_Score.train.'+i+'mer'+'\t'+w_parameter['-Out']+'/Unique_Reads.'+'Combined_Z_Score.train.'+i+'mer'
        process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        function_Run_time(start_time)
        Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
        outfile_list.append(w_parameter['-Out']+'/Unique_Reads.'+'Combined_Z_Score.train.'+i+'mer')

elif(w_parameter['-U'] == '2'):
    folder_path = w_parameter['-Out']+'/Redundant_Reads/'
    for i in nmer_list:
        cmd = 'cp'+'\t'+folder_path+'Combined_Z_Score.train.'+i+'mer'+'\t'+w_parameter['-Out']+'/Redundant_Reads.'+'Combined_Z_Score.train.'+i+'mer'
        process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        function_Run_time(start_time)
        Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
        outfile_list.append(w_parameter['-Out']+'/Redundant_Reads.'+'Combined_Z_Score.train.'+i+'mer')

elif(w_parameter['-U'] == '3'):
    folder_path = w_parameter['-Out']+'/Unique_Reads/'
    for i in nmer_list:
        cmd = 'cp'+'\t'+folder_path+'Combined_Z_Score.train.'+i+'mer'+'\t'+w_parameter['-Out']+'/Unique_Reads.'+'Combined_Z_Score.train.'+i+'mer'
        process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        function_Run_time(start_time)
        Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
        outfile_list.append(w_parameter['-Out']+'/Unique_Reads.'+'Combined_Z_Score.train.'+i+'mer')
    folder_path = w_parameter['-Out']+'/Redundant_Reads/'
    for i in nmer_list:
        cmd = 'cp'+'\t'+folder_path+'Combined_Z_Score.train.'+i+'mer'+'\t'+w_parameter['-Out']+'/Redundant_Reads.'+'Combined_Z_Score.train.'+i+'mer'
        process = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        function_Run_time(start_time)
        Function_Checking_Running_Status(process.returncode, process.stdout + process.stderr)
        outfile_list.append(w_parameter['-Out']+'/Redundant_Reads.'+'Combined_Z_Score.train.'+i+'mer')

print('Training Process finished ...')

print('\n')
print('*** Summary (Training)***')
print('Generated training files: '+str(len(outfile_list)))
print('\n'.join(outfile_list))
print('Those training files will be used for future Aptamer binding potential predictions.')
print('************************')
print('\n')
runningTime = time.time() - start_time
print("Total Running Time:  %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60))