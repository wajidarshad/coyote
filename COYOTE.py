# -*- coding: utf-8 -*-
"""
Created on Sat June 25 11:50:26 2022

@author: Dr. Wajid A. Abbasi
"""
from Bio.Data import IUPACData 
import os,random
import numpy as np
from sklearn import preprocessing



def validate(sequence):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWYXBUZ'
    sequence = sequence.upper()
    #if len(set(sequence).difference(amino_acids)):
       # raise ValueError("Input sequence contains non-standard amino acid codes")
    if len(sequence)<21 or len(set(sequence).difference(amino_acids))>0:
        #raise ValueError("Input sequence is shorter than %d amino acids." %window_size)
        return 0
    else:
        return 1
def standarized_normalized(features, norm='yes'):
    #features=preprocessing.scale(features)
    feats_mean=np.load('features_mean.npy')
    feats_std=np.load('features_std.npy')
    features=(features-feats_mean)/(feats_std+0.0001)
    if norm=='yes':
        features=preprocessing.normalize(features)
    return features
def get_features_proPy_list_Full_des(seq):
        """
        This function takes protein sequence and generate features using
        proPy packege. IN this a search is implemented to get dictionary of all
        terms to get same feature lenght. Normalized fatures will
        be saved at given path with prot id.
        """
        from propy.PyPro import GetProDes

        AA=["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I","L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
        keyList=np.load('key_list_proPy.npy')
        feature=[]
        seq=str(seq)
        seq=seq.replace('X', random.choice(IUPACData.protein_letters))
        seq=seq.replace('B', random.choice(IUPACData.protein_letters))
        seq=seq.replace('Z', random.choice(IUPACData.protein_letters))
        seq=seq.replace('U', random.choice(IUPACData.protein_letters))
        seq=seq.replace('O', random.choice(IUPACData.protein_letters))
        des=GetProDes(seq)
        propy_FD=des.GetALL()
        for key in keyList:
            
            try:
                feature.append(propy_FD[key])
            except KeyError:
                feature.append(0.0)
                #feat=feat/np.linalg.norm(feat)
                #print final_feat_dict,features
                #pdb.set_trace()
        feature=feature/np.linalg.norm(feature)

        return feature

def get_scores(features):
    """
    Obtain raw window scores for all windows
    """
    W_propy=np.load('weight_vector_SVM_propy_glycoproteins_identification.npy')
    bias_propy=np.load('bias_SVM_propy_glycoproteins_identification.npy')
    score =  np.dot(standarized_normalized([features]),W_propy[0])+bias_propy
    return score
def predict_glycoprotein(sequence):
    if not(validate(sequence)):
        print('Input sequence is shorter than 21 amino acids')
        return
    feats=get_features_proPy_list_Full_des(sequence)
    Score=get_scores(feats)
    return Score
