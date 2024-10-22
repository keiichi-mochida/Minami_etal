#!/usr/bin/env python
#coding: UTF-8
#nohup python salmon.py > salmon_out.log & #183

from pathlib import Path
resultsfastqpassname = ["/home/kanatani/20221208BdPDMA/results/"]         #pass***
resultsfastqpass = (" ".join(map(str,resultsfastqpassname)))
shortsamlelist = []

from pathlib import Path
from subprocess import check_call
import os
import re
import pathlib

####salmon_index******************************************
indexnewname = ["Bd21_salmon_index"]                                  #index***
indexname = (" ".join(map(str,indexnewname)))
print("indexname : " + indexname)
#print("{index}".format(index = indexname))

sample_info = "sample_info.txt"                                       #No\t samplename

with open(sample_info) as f:
    next(f)
#    for i in range(1): #test
    for i in range(45):
        lines = f.readline()
        lines_list = lines.split()
        No = lines_list[0]
        check_call("(pigz -d -p 70 {rp}trimmedfastq/trimmed_{S}_1.fq.gz)".format(rp = resultsfastqpass ,S= No ), shell=True)
        check_call("(pigz -d -p 70 {rp}trimmedfastq/trimmed_{S}_2.fq.gz)".format(rp = resultsfastqpass ,S= No ), shell=True)
        check_call("(salmon quant -p 80 -i {rp}/../ref/{index} -l A -1 {rp}trimmedfastq/trimmed_{S}_1.fq -2 {rp}trimmedfastq/trimmed_{S}_2.fq --validateMappings -o {rp}/salmon/{S}.fastq_salmon_quant)".format (index = indexname ,rp = resultsfastqpass ,S= No ),shell=True)
        check_call("(pigz -p 70 {rp}trimmedfastq/trimmed_{S}_1.fq)".format(rp = resultsfastqpass ,S= No ), shell=True)
        check_call("(pigz -p 70 {rp}trimmedfastq/trimmed_{S}_2.fq)".format(rp = resultsfastqpass ,S= No ), shell=True)
