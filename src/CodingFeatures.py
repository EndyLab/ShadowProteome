
class DNAFeatureExtraction:
    '''parent class for dna translation and feature initialization'''
    def __init__(self, sequence_data, list_indices=None):
        self.peptides = None # initialization for subsequent self-checking
        # accepts single 'master' sequence
        if isinstance(sequence_data, str):
            self.seq = sequence_data.upper()
            self.fragments, self.frag_index = stop_delimited_fragments(self.seq)
        # or a pre-split list of fragments.  If no indicies provided, [0]*n used.
        elif isinstance(sequence_data, list):
            self.fragments = [seq.upper() for seq in sequence_data]
            if list_indices:
                self.frag_index = list_indices
            else:    
                self.frag_index = [0]*len(self.fragments)
        
    def decode_dna(self, lookup):
        '''translate dna sequence, compile derived AA features'''
        self.peptides = []
        self.peptide_id = []
        self.codon_pref = []
        self.aa_hydrophob = []
        self.aa_aromatic = []
        self.aa_class = []
        for j,seq in enumerate(self.fragments):
            prot = ''
            pid, pref, aa_h, aa_a, aa_c = [],[],[],[],[]
            for i in range(0,len(seq),3):
                entry = lookup[seq[i:i+3]]
                prot += entry['aa']
                # add'nl properties:
                pid.append(self.frag_index[j])
                pref.append(entry['f_aa']/entry['r_aa'])
                aa_h.append(entry['hydrophob_7'])
                aa_a.append(entry['aromatic'])
                aa_c.append(entry['class'])
            self.peptides.append(prot)
            self.peptide_id.append(pid)
            self.codon_pref.append(pref)
            self.aa_hydrophob.append(aa_h)
            self.aa_aromatic.append(aa_a)
            self.aa_class.append(aa_c)
        return self
    
    def zip_all(self):
        '''list of (amino_acid, codon_pref, hydrophob, aromatic, class, global_peptide_index)'''
        if self.peptides is None:
            print('AttributeError:  Nice try.  You have to decode the DNA first.')
        else:
            return [list(zip(aa,pref,aah,aaa,aac,pix)) 
                    for aa,pref,aah,aaa,aac,pix 
                    in zip(self.peptides,self.codon_pref,self.aa_hydrophob,
                           self.aa_aromatic,self.aa_class,self.peptide_id)]


# DNAFeatureExtraction accessory function
import re
from collections import defaultdict
def stop_delimited_fragments(seq):
    '''Find/Split all possible stop-delimited coding frames\naggregate to single list\nreturns (sequences, global indices)'''
    frames = [0,1,2] # frames
    stops_reg = re.compile('|'.join(['TAA','TAG','TGA'])) # define stops  
    elements = defaultdict(list) # initialize result container
    indicies = defaultdict(list) # initialize indicies container
    idx = {0:0,1:1,2:2} # starting frame indicies

    # progressively search for stop codons, increment result vectors by frame
    i = max(idx.values())
    while i < len(seq):
        x = re.search(stops_reg, seq[i:]) # find next stop
        if x:
            frame = (x.start()+i)%3 # establish frame
            sequence = seq[idx[frame]:x.end()+i] # grab sequence from previous frame end
    #         print(frame, x.group(), sequence)
            elements[frame].append(sequence)
            indicies[frame].append(idx[frame])
            idx[frame] = x.end()+i # update frame start index
            i = max(idx.values())-2 # set new search start (next frame)

        else: # cap all vectors with the remaining non-stopped sequence
            for x in frames:
                sequence = seq[idx[x]:]
                elements[x].append(sequence[:len(sequence)-(len(sequence)%3)]) # trim to full codons
                indicies[x].append(idx[x])
            break
    # aggregate to single list
    fragments = [x for frame in elements.values() for x in frame]
    frag_index = [x for frame in indicies.values() for x in frame] # global index
    return fragments, frag_index


import numpy as np
import pandas as pd
def codon_prefs(seq,w=None):
    '''calculate preference scores for ALL codons'''
    scores = []
    for i,s in enumerate(seq):
        if len(seq[i:i+3])==3:
            entry = lookup[seq[i:i+3]]
            scores.append(entry['f_aa']/entry['r_aa'])
            
    # group by reading frame
    rf_scores = {}
    lmin = 1e6
    for j in [0,1,2]:
        rf_scores[j] = scores[j::3]
        if len(rf_scores[j])<lmin: # find min length
            lmin = len(rf_scores[j])
            
    # truncate to min reading frame length
    for i in rf_scores.keys():
        rf_scores[i] = rf_scores[i][:lmin]
        
    # cast to dataframe
    df = pd.DataFrame.from_dict(rf_scores)
    
    # calculate rolling probability score
    if isinstance(w, int):
        df = df.rolling(window=w,center=True).apply(lambda x: np.exp(sum(np.log(x))/w)) # from the 1983 paper
        
    return df


import pandas as pd
def codon_lookup(path='../data/CodonUsage_ecoli.csv', extrude_=None):
    '''import codon lookup csv as both a pd.DataFrame and lookup dict.  returns (codons,lookup)'''
    codons = pd.read_csv(path)
    codons = codons.set_index('codon')
    if extrude_:
        weights_new = codons.groupby('aa').f_aa.apply(extrude, k=extrude_)
        codons['f_mod'] = [value for probs in weights_new for value in probs]
    # print(codons.head(5))
    lookup = codons.to_dict(orient='index')
    return codons, lookup


# accessory function for codon_lookup.  not generally used.
def extrude(y,k=4):
    x = [(z+.5)**k for z in y]
    xs = [z/sum(x) for z in x]
    return xs
