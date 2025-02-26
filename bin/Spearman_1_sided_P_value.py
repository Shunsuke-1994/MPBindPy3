#!/usr/bin/env python3
import sys
import numpy as np
from scipy.stats import spearmanr

if len(sys.argv) < 3:
    print("Usage: python3 script.py <vector_1> <vector_2>")
    sys.exit(1)

# コマンドライン引数からカンマ区切りの数値列をリストに変換
vector_1 = [float(x) for x in sys.argv[1].split(',')]
vector_2 = [float(x) for x in sys.argv[2].split(',')]

# Spearmanの順位相関検定（二側検定）
rho, p_two = spearmanr(vector_1, vector_2)

# p_twoがNaNの場合は1を返す
if np.isnan(p_two):
    p_value = 1
else:
    # alternative='greater' に合わせる（相関が正の場合は p/2、負の場合は 1 - (p/2)）
    if rho > 0:
        p_value = p_two / 2
    else:
        p_value = 1 - (p_two / 2)

print(p_value)
