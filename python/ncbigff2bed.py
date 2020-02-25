##犹豫ncbi的注释太乱，导致我必须重新生成bed文件与重新填写id

###获取cds与mrna的对应关系

##山羌,ncbi上没有11号染色体

f =open("/home/mwshi/project/xunlu/muntjac/GCA_008787405.1_UCB_Mree_1.0_genomic.gff","r")
chr2id = {"CM018478.1":"1",
          "CM018479.1": "2","CM018480.1":"3","CM018481.1":"4","CM018482.1":"5","CM018483.1":"6",
          "CM018484.1": "7","CM018485.1":"8","CM018486.1":"9","CM018487.1":"10",
          "CM018488.1": "12","CM018489.1":"13","CM018490.1":"14","CM018491.1":"15","CM018492.1":"16",
          "CM018493.1": "17","CM018494.1":"18","CM018495.1":"19","CM018496.1":"20","CM018497.1":"21",
          "CM018498.1": "22","CM018499.1": "23","CM018500.1": "X"}
###
mrna2cds = {}
for line in f:
    if line.startswith("#"):
        continue
    line = line.strip().split("\t")
    if line[2]=="CDS":
        if line[8].startswith("ID"):
            aaaa= line
            cdsid = line[8].split(";")[0]
            mrnaid = line[8].split(";")[1]
            cdsid = cdsid.lstrip("ID=")
            cdsid=cdsid.split(".")[0]
            mrnaid = mrnaid.lstrip("Parent=")
            mrna2cds[mrnaid]=cdsid

###提取bed,id不能带.

f =open("/home/mwshi/project/xunlu/muntjac/GCA_008787405.1_UCB_Mree_1.0_genomic.gff","r")
w = open("/home/mwshi/project/xunlu/muntjac/muntjac.bed","w")
cds2mrna = {}
for line in f:
    if line.startswith("#"):
        continue
    line = line.strip().split("\t")
    if line[2] == "mRNA":
        if line[8].startswith("ID"):
            if line[0] in chr2id.keys():
                chrm = chr2id[line[0]]
            else:
                chrm = line[0]
            start = line[3]
            end = line[4]
            strand = line[6]
            cdsname =  mrna2cds[line[8].split(";")[0].lstrip("ID=rna")]
            w.write(chrm+"\t"+start+"\t"+end+"\t"+cdsname+"\t"+"0"+"\t"+strand+"\n")

w.close()
###fa文件去名
f = open("/home/mwshi/project/xunlu/muntjac/GCA_008787405.1_UCB_Mree_1.0_cds_from_genomic.fna","r")
w = open("/home/mwshi/project/xunlu/muntjac/muntjac.cds","w")
for line in f:
    if line.startswith("#"):
        w.write(line)
    if line.startswith(">"):
        newline = line.strip().split(" ")[0]
        newline = newline.split("_")
        nnn = ">"+newline[1]+"-"+newline[2].split(".")[0]+"\n"
        w.write(nnn)
    else:
        w.write(line)
f.close()
w.close()


##犹豫ncbi的注释太乱，导致我必须重新生成bed文件与重新填写id

###获取cds与mrna的对应关系

##马鹿

f =open("/home/mwshi/project/xunlu/reddeer/GCA_002197005.1_CerEla1.0_genomic.gff","r")
chr2id = {"CM008008.1":"1",
          "CM008009.1": "2","CM008010.1":"3","CM008011.1":"4","CM008012.1":"5","CM008013.1":"6",
          "CM008014.1": "7","CM008015.1":"8","CM008016.1":"9","CM008017.1":"10","CM008018.1":"11",
          "CM008019.1": "12","CM008020.1":"13","CM008021.1":"14","CM008022.1":"15","CM008023.1":"16",
          "CM008024.1": "17","CM008025.1":"18","CM008026.1":"19","CM008027.1":"20","CM008028.1":"21",
          "CM008029.1": "22","CM008030.1": "23","CM008031.1": "24","CM008032.1": "25",
     "CM008033.1": "26","CM008034.1": "27","CM008035.1": "28","CM008036.1": "29","CM008037.1": "30",
    "CM008038.1": "31","CM008039.1": "32","CM008040.1": "33","CM008041.1": "X",
    "CM008042.1": "Y"}
###
mrna2cds = {}
for line in f:
    if line.startswith("#"):
        continue
    line = line.strip().split("\t")
    if line[2]=="CDS":
        if line[8].startswith("ID"):
            aaaa= line
            cdsid = line[8].split(";")[0]
            mrnaid = line[8].split(";")[1]
            cdsid = cdsid.lstrip("ID=")
            cdsid=cdsid.split(".")[0]
            mrnaid = mrnaid.lstrip("Parent=")
            mrna2cds[mrnaid]=cdsid

###提取bed,id不能带.

f =open("/home/mwshi/project/xunlu/reddeer/GCA_002197005.1_CerEla1.0_genomic.gff","r")
w = open("/home/mwshi/project/xunlu/reddeer/reddeer.bed","w")
cds2mrna = {}
for line in f:
    if line.startswith("#"):
        continue
    line = line.strip().split("\t")
    if line[2] == "mRNA":
        if line[8].startswith("ID"):
            if line[0] in chr2id.keys():
                chrm = chr2id[line[0]]
            else:
                chrm = line[0]
            start = line[3]
            end = line[4]
            strand = line[6]
            cdsname =  mrna2cds[line[8].split(";")[0].lstrip("ID=rna")]
            w.write(chrm+"\t"+start+"\t"+end+"\t"+cdsname+"\t"+"0"+"\t"+strand+"\n")

w.close()
###fa文件去名
f = open("/home/mwshi/project/xunlu/reddeer/GCA_002197005.1_CerEla1.0_cds_from_genomic.fna","r")
w = open("/home/mwshi/project/xunlu/reddeer/reddeer.cds","w")
for line in f:
    if line.startswith("#"):
        w.write(line)
    if line.startswith(">"):
        newline = line.strip().split(" ")[0]
        newline = newline.split("_")
        nnn = ">"+newline[1]+"-"+newline[2].split(".")[0]+"\n"
        w.write(nnn)
    else:
        w.write(line)
f.close()
w.close()
