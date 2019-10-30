###读入文件，获取基因名
###shimw
import os
genelist = []
genemap = {}
with open("/home/mwshi/project/xunlu/annotation.txt" , "r") as f:
    for line in f:
        i = line.strip().split("\t")
        id = i[0]
        if i[8]!="-":
            symbol = i[7]
        else:
            symbol = "-"
        genemap[id] = symbol.strip('"')
        genelist.append(symbol)
kkk = []
w = open("/home/mwshi/project/xunlu/id2symbol.txt",'w')
for i in genemap.keys():
    if len(genemap[i].split())==1:
        w.write(i+'\t'+genemap[i]+"\n")

w.close()








