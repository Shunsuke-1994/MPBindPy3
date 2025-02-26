import sys

folder_path = sys.argv[1]  # 例: /datadev/homes/pjiang/work/Aptamer_Pipeline/Peng_Test/Peng_PRAS_Out_U1/Unique_Reads/
RC = sys.argv[2]           # "NULL" またはそれ以外
nmer = sys.argv[3]         # 例: "5,6"

def function_combined_Z_score(Z_score_list):
    Z_score_list = [float(i) for i in Z_score_list]
    k = float(len(Z_score_list))
    combined_Z_score = sum(Z_score_list) / (k ** 0.5)
    return combined_Z_score

def function_combine_Z1_Z2_Z3_Z4(Z1_Z2_file, Z3_file, Z4_file, outfile):
    w_Z1 = {}  # Fisher.Substring
    w_Z2 = {}  # Fisher.Seq
    w_Z3 = {}  # Spearman.Substring
    w_Z4 = {}  # Spearman.Seq
    w_combine_Z = {}
    
    # Z1_Z2 の読み込み
    with open(Z1_Z2_file, 'r') as infile:
        header_temp = infile.readline()
        for line in infile:
            tokens = line.split()
            if not tokens:
                continue
            motif = tokens[0]
            Z1 = tokens[-2]
            Z2 = tokens[-1]
            w_Z1[motif] = Z1
            w_Z2[motif] = Z2
            w_combine_Z[motif] = ''
    
    # Z3 の読み込み
    with open(Z3_file, 'r') as infile:
        header_temp = infile.readline()
        for line in infile:
            tokens = line.split()
            if not tokens:
                continue
            motif = tokens[0]
            Z3 = tokens[-1]
            w_Z3[motif] = Z3
    
    # Z4 の読み込み
    with open(Z4_file, 'r') as infile:
        header_temp = infile.readline()
        for line in infile:
            tokens = line.split()
            if not tokens:
                continue
            motif = tokens[0]
            Z4 = tokens[-1]
            w_Z4[motif] = Z4

    # 各モチーフのCombined Z Scoreを計算し、リストに格納
    all_list = []
    for motif in w_combine_Z.keys():
        Z1 = w_Z1.get(motif, "0")
        Z2 = w_Z2.get(motif, "0")
        Z3 = w_Z3.get(motif, "0")
        Z4 = w_Z4.get(motif, "0")
        combined_Z_this = float(function_combined_Z_score([Z1, Z2, Z3, Z4]))
        out_line = f"{motif}\t{Z1}\t{Z2}\t{Z3}\t{Z4}\t{combined_Z_this}\n"
        all_list.append((combined_Z_this, out_line))
    # 降順にソート
    all_list.sort(key=lambda x: x[0], reverse=True)
    
    with open(outfile, 'w') as outfh:
        header = ('Motif\t'
                  'Z1[Fisher.Substring]\t'
                  'Z2[Fisher.Seq]\t'
                  'Z3[Spearman.Substring]\t'
                  'Z4[Spearman.Seq]\t'
                  'Combined_Z_Score\n')
        outfh.write(header)
        for _, line_str in all_list:
            outfh.write(line_str)

def function_combine_Z3_Z4(Z3_file, Z4_file, outfile):
    w_Z3 = {}  # Spearman.Substring
    w_Z4 = {}  # Spearman.Seq
    w_combine_Z = {}
    
    # Z3 の読み込み
    with open(Z3_file, 'r') as infile:
        header_temp = infile.readline()
        for line in infile:
            tokens = line.split()
            if not tokens:
                continue
            motif = tokens[0]
            Z3 = tokens[-1]
            w_Z3[motif] = Z3
            w_combine_Z[motif] = ''
    
    # Z4 の読み込み
    with open(Z4_file, 'r') as infile:
        header_temp = infile.readline()
        for line in infile:
            tokens = line.split()
            if not tokens:
                continue
            motif = tokens[0]
            Z4 = tokens[-1]
            w_Z4[motif] = Z4

    all_list = []
    for motif in w_combine_Z.keys():
        Z3 = w_Z3.get(motif, "0")
        Z4 = w_Z4.get(motif, "0")
        combined_Z_this = float(function_combined_Z_score([Z3, Z4]))
        out_line = f"{motif}\t{Z3}\t{Z4}\t{combined_Z_this}\n"
        all_list.append((combined_Z_this, out_line))
    all_list.sort(key=lambda x: x[0], reverse=True)
    
    with open(outfile, 'w') as outfh:
        header = ('Motif\t'
                  'Z3[Spearman.Substring]\t'
                  'Z4[Spearman.Seq]\t'
                  'Combined_Z_Score\n')
        outfh.write(header)
        for _, line_str in all_list:
            outfh.write(line_str)

if RC != 'NULL':
    for i in nmer.split(','):
        Z1_Z2_file = folder_path + 'Z1_Z2.' + i + 'mer'
        Z3_file = folder_path + 'Z3.' + i + 'mer'
        Z4_file = folder_path + 'Z4.' + i + 'mer'
        outfile = folder_path + 'Combined_Z_Score.train.' + i + 'mer'
        function_combine_Z1_Z2_Z3_Z4(Z1_Z2_file, Z3_file, Z4_file, outfile)
else:
    for i in nmer.split(','):
        Z3_file = folder_path + 'Z3.' + i + 'mer'
        Z4_file = folder_path + 'Z4.' + i + 'mer'
        outfile = folder_path + 'Combined_Z_Score.train.' + i + 'mer'
        print(Z3_file)
        print(Z4_file)
        function_combine_Z3_Z4(Z3_file, Z4_file, outfile)
