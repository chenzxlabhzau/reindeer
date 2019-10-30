###reindeer mapping
###shimw
###2019/10/17
import os
import subprocess as sub
def qsub(ppn,pbsfile):
    cmd = "qsub -q middle -V -l nodes=1:ppn="+str(ppn)+" "+pbsfile
    sub.call(cmd, shell=True)

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

os.mkdir("/home/mwshi/project/xunlu/bam/pbs")
####建立index
os.mkdir("/home/mwshi/project/xunlu/starindex")
with open("/home/mwshi/project/xunlu/bam/pbs/000_starindex.pbs","w") as w:
    w.write("STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /home/mwshi/project/xunlu/starindex/ --genomeFastaFiles /home/lchen/nmg/02_genome_assembly/xunlu.fa --sjdbGTFfile /home/mwshi/project/genome/reindeer.gtf  --sjdbOverhang 149")

qsub(8,"/home/mwshi/project/xunlu/bam/pbs/000_starindex.pbs")

for i in samples:
    with open("/home/mwshi/project/xunlu/bam/pbs/001_starmapping_"+i+".pbs","w") as w:
        w.write("mkdir -p /home/mwshi/project/prostate/bam/"+i+"\n")
        line = "STAR --runThreadN 10 --twopassMode Basic --readFilesCommand zcat --outReadsUnmapped Fastx --chimSegmentMin 20 --chimJunctionOverhangMin 20 --genomeDir /home/mwshi/project/xunlu/starindex --outSAMtype BAM SortedByCoordinate --readFilesIn "
        line = line + "/home/lchen/nmg/01_rna_raw_data/"+i+"_1.fq.gz "+"/home/lchen/nmg/01_rna_raw_data/" +i+"_2.fq.gz "
        line = line + " --outFileNamePrefix /home/mwshi/project/xunlu/bam/"+i+"/"+i
        w.write(line)


os.mkdir("/home/mwshi/project/xunlu/bam/log")
os.chdir("/home/mwshi/project/xunlu/bam/log")
for i in samples:
    qsub(10,"/home/mwshi/project/xunlu/bam/pbs/001_starmapping_"+i+".pbs")



####count计数

os.mkdir("/home/mwshi/project/xunlu/readcount")
os.mkdir("/home/mwshi/project/xunlu/readcount/pbs")

for i in samples:
    with open("/home/mwshi/project/xunlu/readcount/pbs/002_htcount" + i + ".pbs", "w") as w:
        cmd = "htseq-count -f bam -r pos -s no -t exon -i transcript_id /home/mwshi/project/xunlu/bam/"+i+"/"+i+"Aligned.sortedByCoord.out.bam "
        cmd = cmd + "/home/mwshi/project/genome/reindeer.gtf>/home/mwshi/project/xunlu/readcount/"+i+".txt"
        w.write(cmd)
    qsub(2, "/home/mwshi/project/xunlu/readcount/pbs/002_htcount" + i + ".pbs")

cmd = "Rscript /home/mwshi/project/prostate/script/readcount.R /home/mwshi/project/xunlu/readcount/"
sub.call(cmd, shell=True)

