# -*- coding: utf-8 -*-
"""
Created on Sat June 25 16:19:12 2022

@author: Dr. Wajid A. Abbasi
"""

from Bio import SeqIO
from COYOTE import *

    
    
all_seqs_recs=list(SeqIO.parse("example.fasta", "fasta"))#Give your fasta file instead of example.fasta
print ('Glycoprotein Identification...')
print ('Higher the score more chances of a protein to be a Glycoprotein.')

for rec in all_seqs_recs:
    if validate(rec.seq)==0:
        print('Input sequence is not valid')
    prediction_results=predict_glycoprotein(rec.seq)
    print("The predicted score for the protein "+str(rec.id)+" to be a Glycoprotein is:",prediction_results)

