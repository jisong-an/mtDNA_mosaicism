#!/home/users/anjisong/anaconda3/envs/handygenome/bin/python3.10

################################
#
#   Simulation: random mutation in all regions
#
#   Jisong An
#   jisong0415@kaist.ac.kr
#
################################


## arg1 : MT fasta of the individual
## arg2 : variant summary of the individual (column : context \t n  ;  n : variant count of each context )
## arg3 : output TSV



from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os, sys
import random
from collections import Counter




def make_chunk_3nt(mt_seq, gene_df):

    chunk_merge = []

    for row in gene_df.loc[gene_df[15]=="0,"].itertuples(index=False):
            posstart = int(getattr(row,'_4'))
            posend = int(getattr(row,'_5'))
            strand = str(getattr(row,'_3'))
            if strand == "-":
                x = mt_seq.seq[posstart:posend].reverse_complement()
                chunks = [x[i:i+3] for i in range(0,len(x),3)]
                chunk_merge.extend(chunks)
            else: 
                x = mt_seq.seq[posstart:posend]
                chunks = [x[i:i+3] for i in range(0,len(x),3)]
                chunk_merge.extend(chunks)

    chunks_3nt = [x for x in chunk_merge if len(x)==3]

    return chunks_3nt



def generate_random_mut_make_chunk_3nt(mt_seq, var_df, gene_df):

    mt_seq_temp = mt_seq

    for i,row in var_df.iterrows():

        ref=row.context.split('>')[0]
        alt=row.context.split('>')[1]
        n=int(row.n)  # how many variants in the context

        ref_idx = [i for i in range(len(mt_seq)) if mt_seq.seq.startswith(ref,i)]  # all index of 'ref' base
        random_idx = random.sample(ref_idx,n)  # randomly select index

        for j in random_idx:  # change base into 'alt' at selected index

            base_select = list(mt_seq_temp)
            base_select[j] = alt
            mt_seq_temp = ''.join(base_select)


    mt_seq_temp = Seq(mt_seq_temp)
    chunk_merge_new = []


    # extract coding region
    for row in gene_df.loc[gene_df[15]=="0,"].itertuples(index=False):
            posstart = int(getattr(row,'_4'))
            posend = int(getattr(row,'_5'))
            strand = str(getattr(row,'_3'))
            if strand == "-":
                x = mt_seq_temp[posstart:posend].reverse_complement()
                chunks = [x[i:i+3] for i in range(0,len(x),3)]
                chunk_merge_new.extend(chunks)
            else: 
                x = mt_seq_temp[posstart:posend]
                chunks = [x[i:i+3] for i in range(0,len(x),3)]
                chunk_merge_new.extend(chunks)
    
    # check 'codon' in coding region
    chunk_new_3nt = [x for x in chunk_merge_new if len(x)==3]

    return chunk_new_3nt



def compare_consq(chunk_3nt, chunk_new_3nt):

    total = 0
    syn = 0
    nonsyn = 0
    nonsense = 0 
    stoploss = 0
    syn_var=[]
    nonsyn_var=[]
    nonsense_var=[]
    stoploss_var=[]
    
    # compare consequence
    for x,y in zip(chunk_new_3nt, chunk_3nt):
        if x!=y:
            total+=1
            newaa=x.translate(table=2)  # amino acid
            oldaa=y.translate(table=2)
            
            if oldaa==newaa:
                syn+=1
                syn_var.extend([f'{y[i]}>{x[i]}' for i in range(3) if x[i]!=y[i]])
            elif newaa=="*":
                nonsense+=1
                nonsense_var.extend([f'{y[i]}>{x[i]}' for i in range(3) if x[i]!=y[i]])
            elif oldaa!=newaa:
                nonsyn+=1
                nonsyn_var.extend([f'{y[i]}>{x[i]}' for i in range(3) if x[i]!=y[i]])
            elif oldaa=="*":
                stoploss+=1
                stoploss_var.extend([f'{y[i]}>{x[i]}' for i in range(3) if x[i]!=y[i]])

    # variant context in each consequence
    syn_varlist = '|'.join(f'{key}:{value}' for key,value in Counter(syn_var).items())
    nonsyn_varlist = '|'.join(f'{key}:{value}' for key,value in Counter(nonsyn_var).items())
    nonsense_varlist = '|'.join(f'{key}:{value}' for key,value in Counter(nonsense_var).items())
    stoploss_varlist = '|'.join(f'{key}:{value}' for key,value in Counter(stoploss_var).items())

    return(f'{total};{syn};{nonsyn};{nonsense};{stoploss};{syn_varlist};{nonsyn_varlist};{nonsense_varlist};{stoploss_varlist}')



def main():

    fasta_input = sys.argv[1]
    df_var      = sys.argv[2]
    df_output   = sys.argv[3]

    mt_seq = SeqIO.index(fasta_input, 'fasta')['MT']
    gene_df = pd.read_csv("/home/users/anjisong/reference/MT/refGene_chrMT.txt", sep='\t', header=None)
    var_df = pd.read_csv(df_var, sep='\t')

    # make codon in coding sequence
    chunk_3nt = make_chunk_3nt(mt_seq, gene_df)
    count_dict = {'iter':[],'count':[]}

    # simulation start
    for i in range(0,10000):
        print(f'iteration : {i}')

        # random generation of mutation & extract coding region & make codon
        chunk_new_3nt = generate_random_mut_make_chunk_3nt(mt_seq, var_df, gene_df)
        count_dict['iter'].append(i+1)

        # check consequence of variants
        count_dict['count'].append(compare_consq(chunk_3nt, chunk_new_3nt))

    df = pd.DataFrame(count_dict)
    df[["total","synonymous","nonsynonymous","nonsense","stop_loss","synonymous_var","nonsynonymous_var","nonsense_var","stop_loss_var"]] = df['count'].str.split(";",expand=True)
    df = df.drop("count", axis=1)
    df.to_csv(df_output, index=False, sep='\t')



if __name__ == '__main__':

    main()