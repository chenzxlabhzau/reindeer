
import argparse,os



def main(fastafile,newfastafile):
    f = open(fastafile,"r")
    w = open(newfastafile,"w")
    for line in f:
        if line.startswith("#"):
            w.write(line)
        if line.startswith(">"):
            newline = line.strip().split(" ")[0]
            if newline.startswith(">ENS") and "." in newline:
                newline = newline.split(".")[0]
            w.write(newline+"\n")
        else:
            w.write(line)
    f.close()
    w.close()





if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fastafile", help="the input file you give,must be fasta file", type=str)
    parser.add_argument("newfastafile", help="the output file you give,must be fasta file", type=str)
    args = parser.parse_args()
    fastafile = args.fastafile
    newfastafile = args.newfastafile
    main(fastafile, newfastafile)

