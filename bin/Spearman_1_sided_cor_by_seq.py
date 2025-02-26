#!/usr/bin/env python3
import sys, os, subprocess

# スクリプトの絶対パスを取得
scriptPath = os.path.abspath(os.path.dirname(__file__))

# コマンドライン引数の取得
All_files = sys.argv[1].split(',')
outfile_name = sys.argv[2]

with open(outfile_name, 'w') as outfile:
    header = 'Motif\tFraction_list\tP_value[1-sided]\tZ_Score\n'
    outfile.write(header)
    
    round_list = range(len(All_files))
    
    w_unique_motif = {}
    w_complex = {}
    
    # 各ファイルを読み込み、モチーフとFractionを格納
    for i in round_list:
        with open(All_files[i], 'r') as infile:
            header_temp = infile.readline()
            for line in infile:
                tokens = line.split()
                if not tokens:
                    continue
                motif = tokens[0]
                Fraction = tokens[6]  # *** by_seq ***
                w_unique_motif[motif] = ''
                w_complex[motif + ':' + str(i)] = Fraction
                
    # 各モチーフについてFractionのリストを作成し、コマンドでP値とZスコアを取得
    for k in w_unique_motif.keys():
        outfile.write(k + '\t')
        fraction_list = []
        for i in round_list:
            fraction_list.append(w_complex[k + ':' + str(i)])
        outfile.write(','.join(fraction_list) + '\t')
        
        # 各Fractionはすでに文字列のはずですが、念のため変換
        fraction_list = [str(j) for j in fraction_list]
        round_this_list = [str(j) for j in round_list]
        
        # Spearman_1_sided_P_value.py の実行
        cmd = "python " + os.path.join(scriptPath, "Spearman_1_sided_P_value.py") + " " + ",".join(fraction_list) + " " + ",".join(round_this_list)
        status, output = subprocess.getstatusoutput(cmd)
        P_value = output.strip()
        
        # from_P_value_to_Z_Score_command.py の実行
        cmd = "python " + os.path.join(scriptPath, "from_P_value_to_Z_Score_command.py") + " " + str(P_value)
        status, output = subprocess.getstatusoutput(cmd)
        Z_Score = output.strip()
        
        outfile.write(str(P_value) + '\t' + str(Z_Score) + '\n')
