#!/home/users/anjisong/anaconda3/envs/handygenome/bin/python3.10
################################
#
#   change MT fasta
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
################################

## change MT fasta based on individual's germline variants
## arg 1 : germline variant dataframe
## arg 2 : output file


from Bio import SeqIO
import pandas as pd
import sys,os

fasta_input = "/home/users/anjisong/reference/MT/MT_hg19.fasta"
df_input = sys.argv[1]
fasta_output = sys.argv[2]

df = pd.read_csv(df_input, sep="\t")

with open(fasta_output, "w") as fasta_out:
    mt_seq = SeqIO.index(fasta_input, 'fasta')['MT']

    for row in df[['POS','REF','ALT']].itertuples(index=False):
        pos = int(getattr(row,'POS'))
        alt = str(getattr(row,'ALT'))
        mt_seq.seq = (mt_seq.seq)[:pos-1] + alt + (mt_seq.seq)[pos:]

    SeqIO.write(mt_seq, fasta_out, 'fasta')