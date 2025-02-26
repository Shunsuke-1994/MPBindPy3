#!/usr/bin/env python3
import sys

def function_reverse_complement(seq):
    seq = seq.upper()
    swap_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    newword = ''.join(swap_dict[letter] for letter in reversed(seq))
    return newword

# 引数の取得
Fasta_file = sys.argv[1]
Forward_primer = sys.argv[2].upper()
Reverse_primer = sys.argv[3].upper()
primer_max_mismatch = int(sys.argv[4])
output_file = sys.argv[5]  # 変数名を正しく output_file としています

# Fastaファイルの1件目を読み、Random_Regionの長さを計算
with open(Fasta_file, 'r') as infile:
    line_name = infile.readline()
    line_seq = infile.readline().rstrip()
Random_Region_len = len(line_seq) - len(Forward_primer) - len(Reverse_primer)

####################################
# 2つのシーケンスの長さが等しい前提で、
# 不一致数がmismatch_cutoff以下なら"Yes"、それ以外は"No"を返す
def function_less_or_eq_mismatch_logical(seq_1, seq_2, mismatch_cutoff):
    if len(seq_1) != len(seq_2):
        print('Error... Len_1 != Len2')
    N_mismatch = 0
    seq_len = len(seq_1)
    for i in range(seq_len):
        if N_mismatch <= mismatch_cutoff:
            if seq_1[i] != seq_2[i]:
                N_mismatch += 1
        else:
            return "No"
    if N_mismatch > mismatch_cutoff:
        return "No"
    return "Yes"

with open(Fasta_file, 'r') as infile, open(output_file, 'w') as outfile:
    # Fasta形式なので、2行ごとに1件のレコードとして処理
    while True:
        line_name = infile.readline()
        if not line_name:
            break
        line_seq = infile.readline()
        if not line_seq:
            break
        seq_whole = line_seq.rstrip().upper()
        # 不明な塩基やドットが含まれていない場合
        if 'N' not in seq_whole and '.' not in seq_whole:
            # フォワードプライマーでチェック
            primer_region = seq_whole[:len(Forward_primer)]
            random_region = seq_whole[len(Forward_primer):len(Forward_primer) + Random_Region_len]
            if function_less_or_eq_mismatch_logical(primer_region, Forward_primer, primer_max_mismatch) == 'Yes':
                outfile.write(random_region.upper() + '\n')
            else:
                # フォワードプライマーが不適合ならリバースプライマーでチェック
                primer_region = seq_whole[:len(Reverse_primer)]
                random_region = seq_whole[len(Reverse_primer):len(Reverse_primer) + Random_Region_len]
                if function_less_or_eq_mismatch_logical(primer_region, Reverse_primer, primer_max_mismatch) == 'Yes':
                    outfile.write(function_reverse_complement(random_region.upper()) + '\n')
        # 該当しない場合は何もしない
