#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 21 11:33:57 2021

@author: rafael
"""
# BIBLIOTECAS
import Bio
# import heapq
import pylab
import urllib
import pandas as pd
import nglview as nv

from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
#from Bio.Alphabet import IUPAC
from collections import Counter
from Bio.Data import CodonTable
from Bio import SeqIO, SearchIO
from Bio.PDB import PDBParser,MMCIFParser
from Bio.SeqUtils import GC,molecular_weight
# from dna_features_viewer import GraphicFeature, GraphicRecord
#from Bio.Alphabet import generic_dna,generic_rna,generic_protein

#IMPORTAÇÃO DA SEQUÊNCIA
#sequencia = SeqIO.read("sequence.fasta", "fasta")

for record in SeqIO.parse("sequence.fasta","fasta"):
    print(record)

#EXPLORAÇÃO DA SEQUÊNCIA
sequencia_do_arquivo = record.seq
sequencia_do_arquivo
len(sequencia_do_arquivo)
molecular_weight(sequencia_do_arquivo)
sequencia_do_arquivo.find('AGA')
sequencia_do_arquivo[0:10].complement()
sequencia_do_arquivo[0:10].reverse_complement()

#MANIPULAÇÃO DA SEQUÊNCIA

sequencia_do_arquivo[0:3]
sequencia_do_arquivo[0:3] + sequencia_do_arquivo[-3:]
sequencia_do_arquivo.find('AGA')
(sequencia_do_arquivo.count('G') + sequencia_do_arquivo.count('C'))/(len(sequencia_do_arquivo)) * 100
sequencia_do_arquivo[0:6]
sequencia_do_arquivo[0:6].complement()
sequencia_do_arquivo[0:6].reverse_complement()

#TRANSCRIÇÃO E TRADUÇÃO

mRNA = sequencia_do_arquivo.transcribe()
mRNA[:10]
mRNA.back_transcribe()[:10]
print(CodonTable.unambiguous_dna_by_id[1])
protein_seq = sequencia_do_arquivo.translate()
protein_seq[:10]
len(protein_seq)
len(sequencia_do_arquivo)
protein_seq[-10:]
protein_seq.back_transcribe()[-10:]
print(CodonTable.unambiguous_dna_by_id[1])
common_amino = Counter(protein_seq)
common_amino.most_common(10)
del common_amino['*']
pylab.bar(common_amino.keys(),common_amino.values())
pylab.title("%i protein sequences\n frquency %i to %i" 
            % (len(common_amino.values()), 
               min(common_amino.values()), 
               max(common_amino.values())))
pylab.xlabel("Amino acid")
pylab.ylabel("frequency")
pylab.show()
protein_list = [str(i) for i in protein_seq.split('*')]
protein_list[:10]
large_proteins = [x for x in protein_list if len(x)> 10]
df = pd.DataFrame({'protein_seq':large_proteins})
df['length'] = df['protein_seq'].apply(len)
df.head()
df.sort_values(by = ['length'], ascending = False)[:10]
one_large_protein = df.nlargest(1,'length')
single_prot = one_large_protein.iloc[0,0]
single_prot
with open("single_seq.fasta","w") as file:
    file.write(">largest_seq \n"+single_prot)

#ALINHAMENTO LOCAL COM NCBI-BLAST

read = SeqIO.read("single_seq.fasta", "fasta")
read.seq
#%%time
result_handle = NCBIWWW.qblast("blastp","pdb",read.seq)
blast_qresult = SearchIO.read(result_handle, "blast-xml")
print(blast_qresult[0:5])
seqid = blast_qresult[0]
details = seqid[0]
print(f"\
Sequence ID:{seqid.id}\n\
description:{seqid.description}\n\
E value:    {details.evalue} \n\
Bit Score:  {details.bitscore}\n\
")
print(f"alignment:\n{details.aln}")

#LENDO ARQUIVO PDB

seqid.id
seqid.id.split('|')[1]
urllib.request.urlretrieve('https://files.rcsb.org/download/6YYT.pdb',
                           '6YYT.pdb')
parser = PDBParser()
structure = parser.get_structure("6YYT","6YYT.pdb")
structure
for chain in structure[0]:
    print(f"chainid: {chain.id}")

#VISUALIZANDO ESTRUTURA DA PROTEINA
nv.demo()
view = nv.show_biopython(structure)
view
view.render_image()
view = nv.show_biopython(structure, gui=True)
view
view.render_image()


