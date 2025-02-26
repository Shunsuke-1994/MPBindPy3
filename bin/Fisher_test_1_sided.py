#!/usr/bin/env python3
import sys
from scipy.stats import fisher_exact

if len(sys.argv) < 5:
    print("Usage: {} a1 a2 b1 b2".format(sys.argv[0]))
    sys.exit(1)

# コマンドライン引数の取得と整数変換
a1 = int(sys.argv[1])
a2 = int(sys.argv[2])
b1 = int(sys.argv[3])
b2 = int(sys.argv[4])

# Rの cbind(c(a1,a2),c(b1,b2)) に相当する2×2分割表を作成
table = [[a1, b1],
         [a2, b2]]

# 1-sided (greater) のFisherの正確検定を実施
_, p_value = fisher_exact(table, alternative='greater')

print(p_value)
