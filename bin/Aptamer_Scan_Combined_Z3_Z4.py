import sys

Train_file = sys.argv[1]
Aptamer_file = sys.argv[2]
Sort_logical = sys.argv[3]  # 'TRUE' or 'FALSE'
Output_file = sys.argv[4]

def function_combined_Z_score(Z_score_list):
    Z_score_list = [float(i) for i in Z_score_list]
    k = float(len(Z_score_list))
    combined_Z_score = sum(Z_score_list) / (k ** 0.5)
    return combined_Z_score

def function_aptamer_scan_meta_Z_score(w_motif_Z_score, aptamer_seq):
    # Python 3ではdict.keys()はビューを返すため、リスト化して先頭のキーを取得
    motif = list(w_motif_Z_score.keys())[0]
    motif_len = len(motif)
    substring_Z_score_list = []
    for i in range(len(aptamer_seq) - motif_len + 1):
        substring = aptamer_seq[i:i+motif_len]
        substring_Z_score_list.append(float(w_motif_Z_score[substring]))
    meta_Z_score = function_combined_Z_score(substring_Z_score_list)
    substring_Z_score_list = [str(round(i, 2)) for i in substring_Z_score_list]
    return ','.join(substring_Z_score_list) + '\t' + str(round(meta_Z_score, 2))

# モチーフZスコアの読み込み
w_motif_Z3 = {}
w_motif_Z4 = {}
w_motif_Z_Combined = {}

with open(Train_file, 'r') as infile:
    header_temp = infile.readline()  # ヘッダ行を読み飛ばす
    for line in infile:
        fields = line.split()
        motif = fields[0]
        Z3 = fields[1]
        Z4 = fields[2]
        Combined_Z_Score = fields[3]
        w_motif_Z3[motif] = Z3
        w_motif_Z4[motif] = Z4
        w_motif_Z_Combined[motif] = Combined_Z_Score

if Sort_logical.upper() == 'FALSE':
    with open(Aptamer_file, 'r') as infile, open(Output_file, 'w') as outfile:
        header = ('Aptamer.Seq\t'
                  'Z3.Scan\tZ3.MetaScore\t'
                  'Z4.Scan\tZ4.MetaScore\t'
                  'Z_Combined.Scan\tZ_Combined.MetaScore\n')
        outfile.write(header)
        for line in infile:
            Seq = line.split()[0].upper()
            outfile.write(Seq + '\t')
            outfile.write(function_aptamer_scan_meta_Z_score(w_motif_Z3, Seq) + '\t')
            outfile.write(function_aptamer_scan_meta_Z_score(w_motif_Z4, Seq) + '\t')
            outfile.write(function_aptamer_scan_meta_Z_score(w_motif_Z_Combined, Seq) + '\n')
elif Sort_logical.upper() == 'TRUE':
    out_list = []
    with open(Aptamer_file, 'r') as infile:
        for line in infile:
            Seq = line.split()[0].upper()
            this_Z3 = function_aptamer_scan_meta_Z_score(w_motif_Z3, Seq)
            this_Z4 = function_aptamer_scan_meta_Z_score(w_motif_Z4, Seq)
            this_combined = function_aptamer_scan_meta_Z_score(w_motif_Z_Combined, Seq)
            # スコアはタブ区切りの最後の要素
            this_score = function_aptamer_scan_meta_Z_score(w_motif_Z_Combined, Seq).split()[-1]
            out_list.append([float(this_score),
                             Seq + '\t' + this_Z3 + '\t' + this_Z4 + '\t' + this_combined + '\n'])
    # 降順ソート
    out_list.sort(reverse=True)
    with open(Output_file, 'w') as outfile:
        header = ('Aptamer.Seq\t'
                  'Z3.Scan\tZ3.MetaScore\t'
                  'Z4.Scan\tZ4.MetaScore\t'
                  'Z_Combined.Scan\tZ_Combined.MetaScore\n')
        outfile.write(header)
        for _, line_str in out_list:
            outfile.write(line_str)
else:
    print('Error ...')
    print('Check -sort option ...')
    sys.exit(1)
