###读入文件，获取基因名
###shimw
import os
genelist = []
with open("/home/mwshi/project/xunlu/OR.txt" , "r") as f:
    for line in f:
        i = line.strip().split("\t")
        # if i[5].split(" ")[5].split("=")[0]=="GN":
        #     gene_name = i[5].split(" ")[5].split("=")[1]
        # else:
        #     print(i)
        for item in i:
            if i[5].startswith("Olfactory"):





