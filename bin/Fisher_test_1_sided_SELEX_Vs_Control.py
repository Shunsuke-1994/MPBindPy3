import sys
import os
import subprocess

# スクリプトの絶対パスを取得
scriptPath = os.path.abspath(os.path.dirname(__file__))

# コマンドライン引数
SELEX_file = sys.argv[1]
Control_file = sys.argv[2]
Outfile_name = sys.argv[3]

###################################################
def function_SELEX_Vs_Control_Fisher_1_sided(SELEX_file, Control_file, Outfile_name):
    w_motif_SELEX = {}  
    w_motif_Control = {}

    # SELEXファイルの読み込み
    with open(SELEX_file, 'r') as infile:
        header_temp = infile.readline()  # ヘッダー行
        for line in infile:
            line = line.strip()
            if not line:
                continue
            tokens = line.split()
            Motif = tokens[0]
            Substring_counts = tokens[1]
            Substring_others = str(int(tokens[2]) - int(tokens[1]))
            Substring_fraction = tokens[3]
            Seq_counts = tokens[4]
            Seq_others = str(int(tokens[5]) - int(tokens[4]))
            Seq_fraction = tokens[6]
            w_motif_SELEX[Motif] = [Substring_counts, Substring_others, Substring_fraction, Seq_counts, Seq_others, Seq_fraction]

    # Controlファイルの読み込み
    with open(Control_file, 'r') as infile:
        header_temp = infile.readline()  # ヘッダー行
        for line in infile:
            line = line.strip()
            if not line:
                continue
            tokens = line.split()
            Motif = tokens[0]
            Substring_counts = tokens[1]
            Substring_others = str(int(tokens[2]) - int(tokens[1]))
            Substring_fraction = tokens[3]
            Seq_counts = tokens[4]
            Seq_others = str(int(tokens[5]) - int(tokens[4]))
            Seq_fraction = tokens[6]
            w_motif_Control[Motif] = [Substring_counts, Substring_others, Substring_fraction, Seq_counts, Seq_others, Seq_fraction]

    # 出力ファイルに結果を書き出し
    with open(Outfile_name, 'w') as outfile:
        header = ('Motif\t'
                  'P_value.Substring[1-sided]\t'
                  'P_value.Seq[1-sided]\t'
                  'Z-Score.Substring\t'
                  'Z-Score.Seq\n')
        outfile.write(header)
        
        for k in w_motif_SELEX.keys():
            outfile.write(k + '\t')
            # Substringの値取得
            Substring_a1 = w_motif_SELEX[k][0]
            Substring_a2 = w_motif_SELEX[k][1]
            Substring_b1 = w_motif_Control[k][0]
            Substring_b2 = w_motif_Control[k][1]
            # Rスクリプト実行（Fisherの1-sidedテスト）
            cmd = ' '.join([
                "python",
                os.path.join(scriptPath, "Fisher_test_1_sided.py"),
                str(Substring_a1),
                str(Substring_a2),
                str(Substring_b1),
                str(Substring_b2)
            ])
            proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            Substring_P_value = proc.stdout.strip()
            
            # P値からZスコアへの変換
            cmd = ' '.join([
                "python",
                os.path.join(scriptPath, "from_P_value_to_Z_Score_command.py"),
                str(Substring_P_value)
            ])
            proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            Substring_Z_score = proc.stdout.strip()        

            # Seqの値取得
            Seq_a1 = w_motif_SELEX[k][3]
            Seq_a2 = w_motif_SELEX[k][4]
            Seq_b1 = w_motif_Control[k][3]
            Seq_b2 = w_motif_Control[k][4]
            # Rスクリプト実行（Fisherの1-sidedテスト）
            cmd = ' '.join([
                "python",
                os.path.join(scriptPath, "Fisher_test_1_sided.py"),
                str(Seq_a1),
                str(Seq_a2),
                str(Seq_b1),
                str(Seq_b2)
            ])
            proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            Seq_P_value = proc.stdout.strip()
            
            # P値からZスコアへの変換
            cmd = ' '.join([
                "python",
                os.path.join(scriptPath, "from_P_value_to_Z_Score_command.py"),
                str(Seq_P_value)
            ])
            proc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            Seq_Z_score = proc.stdout.strip()

            outfile.write(f"{Substring_P_value}\t{Seq_P_value}\t{Substring_Z_score}\t{Seq_Z_score}\n")

#############################################################################################
function_SELEX_Vs_Control_Fisher_1_sided(SELEX_file, Control_file, Outfile_name)
