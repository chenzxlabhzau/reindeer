###由于ncbi自己的注释太烂，导致我只能用先从gff拿到基因对应的蛋白id
### 然后从pep文件中拿到所有的蛋白长度，然后提取所有基因的最长的蛋白序列
from Bio import SeqIO
import os


def get_key (dict, value):
    return [k for k, v in dict.items() if v == value]

def writenewfile(gene2cds,pepfile,outfile):
    ##开始做蛋白长度表
    cdslendict ={}
    for record in SeqIO.parse(pepfile, "fasta"):
        cdslendict[record.id]=len(str(record.seq))
    ##取最长
    gene2maxLcds = {}
    for i in gene2cds.keys():
        maxlencds=""
        genelen = 0
        if len(gene2cds[i])==1 and (gene2cds[i][0] not in cdslendict.keys()):
            continue
        ###跳过ENSBTAG00000024549
        if not all([word in cdslendict.keys() for word in gene2cds[i]]):
            continue
        for j in gene2cds[i]:
            if cdslendict[j]>genelen:
                genelen = cdslendict[j]
                maxlencds = j
        gene2maxLcds[i]=maxlencds
    ###然后开始写pep文件
    with open(outfile,"w") as f:
        for record in SeqIO.parse(pepfile, "fasta"):
            if record.id in gene2maxLcds.values():
                f.write(">"+str(record.id)+"\n")
                f.write(str(record.seq) + "\n")


def ncbimodel(gff,pepfile,outfile):
    #gff="/home/mwshi/project/xunlu/mcl/genome/reindeer/reindeer.gff"
    #pepfile="/home/mwshi/project/xunlu/mcl/genome/reindeer/reindeer.pep"
    #outfile="/home/mwshi/project/xunlu/mcl/genome/reindeer/reindeer.L.pep.fa"
    gfffile = open(gff,"r")

###获取对应关系
    genedict = {}
    cds2mrna = {}
    for line in gfffile:
        if line.startswith("#"):
            continue
        line = line.strip().split("\t")
        if line[2] == "mRNA":
            if line[8].startswith("ID"):
                mrnaid = line[8].split(";")[0]
                geneid = line[8].split(";")[1]
                mrnaid = mrnaid.replace("ID=",'')
                geneid=geneid.replace("Parent=",'')
                genedict[mrnaid]=geneid
        if line[2] == "CDS":
            if line[8].startswith("ID"):
                cdsid = line[8].split(";")[0]
                mrnaid = line[8].split(";")[1]
                cdsid = cdsid.replace("ID=", '')
                ###判断是否是cds-xxxx格式
                if "-" in cdsid:
                    cdsid=cdsid.split("-")[1]
                mrnaid = mrnaid.replace("Parent=", '')
                cds2mrna[cdsid]=mrnaid

    ###z这里遇到了特殊情况，就是有基因没有mrna，但是他有cds，导致没有cds对应的mrna这种情况下，其对应的mrna本身就是其gene
    ###所以这里直接把mrna对到基因上
    ###遇到了mrna的行部位mRNA而是C_gene_segment或者V_gene_segment等，发现他们都只有一个转录本，所以直接跳过
    cds2gene = cds2mrna
    for i in cds2mrna.keys():
        if cds2mrna[i] in genedict.keys():
            cds2gene[i]=genedict[cds2mrna[i]]
    ##然后反向对应id
    gene2cds={}
    for i in set(cds2gene.values()):
        gene2cds[i]=get_key(cds2gene,i)
    writenewfile(gene2cds, pepfile, outfile)


###EVM注释格式
def EVMmode(gff,pepfile,outfile):
    #gff="/home/mwshi/project/xunlu/mcl/genome/reindeer/reindeer.gff"
    #pepfile="/home/mwshi/project/xunlu/mcl/genome/reindeer/reindeer.pep"
    #outfile="/home/mwshi/project/xunlu/mcl/genome/reindeer/reindeer.L.pep.fa"
    gfffile = open(gff, "r")
    ###驯鹿使用EVM注释，mrna中直接就是蛋白id，或者本身没有转录本和蛋白id的区别。获取对应关系
    cds2gene = {}
    for line in gfffile:
        if line.startswith("#"):
            continue
        line = line.strip().split("\t")
        if line[2] == "mRNA":
            if line[8].startswith("ID"):
                cds = line[8].split(";")[0]
                geneid = line[8].split(";")[1]
                cds = cds.replace("ID=",'')
                geneid=geneid.replace("Parent=",'')
                cds2gene[cds]=geneid
    ###翻转id对应关系
    gene2cds={}
    for i in set(cds2gene.values()):
        gene2cds[i]=get_key(cds2gene,i)
    writenewfile(gene2cds, pepfile, outfile)


def ncbiNORmodel(gff,pepfile,outfile,pre):
    #gff="/home/mwshi/project/xunlu/mcl/genome/reddeer/GCA_002197005.1_CerEla1.0_genomic.gff"
    #pepfile="/home/mwshi/project/xunlu/mcl/genome/reddeer/reddeer.pep.fa"
    #outfile="/home/mwshi/project/xunlu/mcl/genome/reddeer/reddeer.L.pep.fa"
    #pre = "cds-"
    gfffile = open(gff,"r")

###获取对应关系
    genedict = {}
    cds2mrna = {}
    for line in gfffile:
        if line.startswith("#"):
            continue
        line = line.strip().split("\t")
        if line[2] == "mRNA":
            if line[8].startswith("ID"):
                mrnaid = line[8].split(";")[0]
                geneid = line[8].split(";")[1]
                mrnaid = mrnaid.replace("ID=",'')
                geneid=geneid.replace("Parent=",'')
                genedict[mrnaid]=geneid
        if line[2] == "CDS":
            if line[8].startswith("ID"):
                cdsid = line[8].split(";")[0]
                mrnaid = line[8].split(";")[1]
                cdsid = cdsid.replace("ID=", '')
                ###判断是否是cds-xxxx格式,去掉pre
                cdsid=cdsid.replace(pre, '')
                mrnaid = mrnaid.replace("Parent=", '')
                cds2mrna[cdsid]=mrnaid
    gfffile.close()
    ###z这里遇到了特殊情况，就是有基因没有mrna，但是他有cds，导致没有cds对应的mrna这种情况下，其对应的mrna本身就是其gene
    ###所以这里直接把mrna对到基因上
    ###遇到了mrna的行部位mRNA而是C_gene_segment或者V_gene_segment等，发现他们都只有一个转录本，所以直接跳过
    cds2gene = cds2mrna
    for i in cds2mrna.keys():
        if cds2mrna[i] in genedict.keys():
            cds2gene[i]=genedict[cds2mrna[i]]
    ##然后反向对应id
    gene2cds={}
    for i in set(cds2gene.values()):
        gene2cds[i]=get_key(cds2gene,i)
    writenewfile(gene2cds, pepfile, outfile)

###WhitetailedDeer
ncbiNORmodel("/home/mwshi/project/xunlu/mcl/genome/WhitetailedDeer/GCF_002102435.1_Ovir.te_1.0_genomic.gff",
             "/home/mwshi/project/xunlu/mcl/genome/WhitetailedDeer/WhitetailedDeer.pep.fa",
             "/home/mwshi/project/xunlu/mcl/genome/WhitetailedDeer/WhitetailedDeer.pepL.fa",
             "cds-")
###reindeer
EVMmode("/home/mwshi/project/xunlu/mcl/genome/reindeer/reindeer.gff",
        "/home/mwshi/project/xunlu/mcl/genome/reindeer/reindeer.pep.fa",
        "/home/mwshi/project/xunlu/mcl/genome/reindeer/reindeer.pepL.fa")

##reddeer
ncbiNORmodel("/home/mwshi/project/xunlu/mcl/genome/reddeer/GCA_002197005.1_CerEla1.0_genomic.gff",
             "/home/mwshi/project/xunlu/mcl/genome/reddeer/reddeer.pep.fa",
             "/home/mwshi/project/xunlu/mcl/genome/reddeer/reddeer.pepL.fa",
             "cds-")
##pig
ncbiNORmodel("/home/mwshi/project/xunlu/mcl/genome/pig/Sus_scrofa.Sscrofa11.1.98.gff3",
             "/home/mwshi/project/xunlu/mcl/genome/pig/pig.pep.fa",
             "/home/mwshi/project/xunlu/mcl/genome/pig/pig.pepL.fa",
             "CDS:")
##muntjac
ncbiNORmodel("/home/mwshi/project/xunlu/mcl/genome/muntjac/GCA_008787405.1_UCB_Mree_1.0_genomic.gff",
             "/home/mwshi/project/xunlu/mcl/genome/muntjac/muntjac.pep.fa",
             "/home/mwshi/project/xunlu/mcl/genome/muntjac/muntjac.pepL.fa",
             "cds-")
###milu,无基因id，所以无法提取基因最长转录本
os.popen('cp /home/mwshi/project/xunlu/mcl/genome/Milu/Milu.pep.fa /home/mwshi/project/xunlu/mcl/genome/Milu/Milu.pepL.fa')

##human
ncbiNORmodel("/home/mwshi/project/xunlu/mcl/genome/human/Homo_sapiens.GRCh38.98.gff3",
             "/home/mwshi/project/xunlu/mcl/genome/human/human.pep.fa",
             "/home/mwshi/project/xunlu/mcl/genome/human/human.pepL.fa",
             "CDS:")
##goat
ncbiNORmodel("/home/mwshi/project/xunlu/mcl/genome/goat/Capra_hircus.ARS1.98.gff3",
             "/home/mwshi/project/xunlu/mcl/genome/goat/goat.pep.fa",
             "/home/mwshi/project/xunlu/mcl/genome/goat/goat.pepL.fa",
             "CDS:")
##giraffe,无基因只有mrna
os.popen('cp /home/mwshi/project/xunlu/mcl/genome/giraffe/giraffe.pep.fa /home/mwshi/project/xunlu/mcl/genome/giraffe/giraffe.pepL.fa')

##Forestmuskdeer,无基因只有mrna
os.popen('cp /home/mwshi/project/xunlu/mcl/genome/Forestmuskdeer/Forestmuskdeer.pep.fa /home/mwshi/project/xunlu/mcl/genome/Forestmuskdeer/Forestmuskdeer.pepL.fa')

##cow.牛的gff中的mrna在pep文件中不存在，导致必须修改函数进行判断是否存在
ncbiNORmodel("/home/mwshi/project/xunlu/mcl/genome/cow/Bos_taurus.ARS-UCD1.2.98.gff3",
             "/home/mwshi/project/xunlu/mcl/genome/cow/cow.pep.fa",
             "/home/mwshi/project/xunlu/mcl/genome/cow/cow.pepL.fa",
             "CDS:")

