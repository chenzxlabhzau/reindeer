###circadian 序列比对
###前期的基因挑选在电脑上的R进行的
import argparse
import os
from Bio import SeqIO
import pandas as pd
import subprocess as sub
def qsub(ppn,pbsfile):
    cmd = "qsub -q middle -V -l nodes=1:ppn="+str(ppn)+" "+pbsfile
    sub.call(cmd, shell=True)


def align_circadian(dir,fastafile,species):
    cmd = "python /home/mwshi/github/reindeer/python/fastaRename.py "+fastafile+ " "+dir+species+".cds"
    sub.call(cmd, shell=True)
    species_circadian = pd.read_csv(dir+species+"_Circadian.csv")

    with open(dir+species+"_Circadian.cds","w") as f:
        for record in SeqIO.parse(dir+species+".cds", "fasta"):
            if record.id in species_circadian["ensembl_transcript_id"].values:
                symbol = species_circadian[species_circadian["ensembl_transcript_id"] == str(record.id)]["external_gene_name"].values[0]
                f.write(">"+str(record.id)+"|"+symbol+"\n")
                f.write(str(record.seq) + "\n")

    os.chdir('/home/mwshi/project/xunlu/human/')

    cmd = "lastal -E 0.05 human_CircadianDB "+dir+species+"_Circadian.cds > "+dir+"human."+species+".maf"
    sub.call(cmd, shell=True)


align_circadian("/home/mwshi/project/xunlu/mouse/","/home/mwshi/project/xunlu/mouse/Mus_musculus.GRCm38.cds.all.fa","mouse")
align_circadian("/home/mwshi/project/xunlu/goat/","/home/mwshi/project/xunlu/goat/Capra_hircus.ARS1.cds.all.fa","goat")
align_circadian("/home/mwshi/project/xunlu/cow/","/home/mwshi/project/xunlu/cow/Bos_taurus.ARS-UCD1.2.cds.all.fa","cow")

###以上常熟失败，使用muscle比对
###重新获取cds序列
def getsquence(species,datapath):
    dataframe = pd.read_csv(datapath)
    seqdict = {}
    with open("/home/mwshi/project/xunlu/"+species+"/"+species+"_L_Circadian.cds","w") as f:
        for record in SeqIO.parse("/home/mwshi/project/xunlu/"+species+"/"+species+".cds", "fasta"):
            if record.id in dataframe[dataframe.columns[0]].values:
                symbol = dataframe[dataframe[dataframe.columns[0]] == str(record.id)][dataframe.columns[1]].values[0]
                f.write(">"+str(record.id)+"|"+symbol+"\n")
                f.write(str(record.seq) + "\n")
                seqdict[str(record.id)] = str(record.seq)
    return seqdict

huamnseq = getsquence("human","/home/mwshi/project/xunlu/MUSCLE/human_L_Circadian.csv")

for i in dataframe["ensembl_transcript_id"]:
    if i not in huamnseq.keys():
        print(i)

