#!/usr/bin/env python3
import sys
from scipy.stats import norm

def function_z_score(p_value):
    if p_value < 10**(-15):
        z_score = 8.0
    elif p_value > (1 - 10**(-15)):
        z_score = -8.0
    else:
        z_score = norm.ppf(1 - p_value)
    return z_score

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python3 from_P_value_to_Z_Score_example.py <P-value>")
        sys.exit(1)
    p_value = float(sys.argv[1])
    z_score = function_z_score(p_value)
    print(z_score)
