#!/usr/bin/env python
import sys
umis = [int(x) for x in sys.stdin.read().split(",")]

umis = sorted(umis)
print(umis[int(len(umis) * 1)-1])  # 40th percentile
