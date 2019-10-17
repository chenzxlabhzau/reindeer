###reindeer mapping
###shimw
###2019/10/17
import os
samples = []
for i in os.listdir("/home/lchen/nmg/01_rna_raw_data/"):
    i = i.split("_")
    if i[0]=="subcutaneous":
        sample = i[0]+"_"+i[1]
    else:
        sample = i[0]
    samples.append(sample)

    samples = list(set(samples))

for i in samples:
    os.mkdir("/home/mwshi/project/xunlu/bam/"+i)




