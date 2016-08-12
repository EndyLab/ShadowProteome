
class DNAFeatureExtraction:
    '''parent class for dna translation and feature initialization'''
    def __init__(self, sequence_data, chunk=True, list_indices=None):
        self.peptides = None # initialization for subsequent self-checking
        # accepts single 'master' sequence
        if isinstance(sequence_data, str):
            self.seq = sequence_data.upper()
            self.fragments, self.frag_index = fragment_stop_delimited(self.seq)
        # or a pre-split list of fragments.  If no indicies provided, [0]*n used.
        elif isinstance(sequence_data, list):
            self.fragments = [seq.upper() for seq in sequence_data]
            if list_indices:
                self.frag_index = list_indices
            else:    
                self.frag_index = [0]*len(self.fragments)
        if chunk:
            self.fragments, self.frag_index = fragment_windowed(zip(self.fragments, self.frag_index))
            

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
def fragment_stop_delimited(seq):
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
    

# DNAFeatureExtraction accessory function
def fragment_windowed(seqs_idx, window=150, offset=60):
    '''takes list of tuple(aa_sequence,index) and returns lists of windows & updated indicies'''
    chunked = [(frag[idn:idn+window], idg+idn) for frag,idg in seqs_idx for idn in range((len(frag)-window)%offset, len(frag)-window+1, offset)]
    windows, wind_index = list(zip(*chunked))
    return list(windows), list(wind_index)


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


## ACCESSORY STATPACK FUCTIONS
def zerocross_pack(X): # intended for hydrophob array
    '''Takes 1d array of numerical values and returns a (5,) np.array of distribution features'''
    cross = []
    for k in range(1,len(X)):
        if X[k]*X[k-1]<0:
            cross.append(k)
    t = np.ediff1d(cross)
    tm = t.mean()
    tsd = t.std()
    tp1 = np.percentile(t,5)
    tp2 = np.percentile(t,95)
    return np.hstack([len(cross)/len(X),tm,tsd,tp1,tp2])

    
def codon_bias_rms(X, lookup, select='FULL'):
    '''Takes a DNA sequence and calculates the RMS[or 0] of codon bias for ALL(1,), [SELECT](1,), or FULL(21,) amino acids'''
    AAs = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    x = [(lookup[X[i:i+3]]['aa'],lookup[X[i:i+3]]['f_aa']-lookup[X[i:i+3]]['r_aa']) for i in range(0,len(X),3)]
    
    if select=='ALL':
        return (sum([z[1]**2 for z in x])/len(x))**.5
    elif select=='FULL':
        df = pd.DataFrame(x)
        values = df[1].groupby(df[0]).apply(lambda x: (sum(x**2)/len(x))**.5)
        return np.hstack([(sum([z[1]**2 for z in test])/len(test))**.5, values[AAs].fillna(0)])
    else:
        tmp = [z[1]**2 for z in x if z[0]==select]
        return (sum(tmp)/len(tmp)) if len(tmp)>0 else 0.


def aa_pack(X, AAs): # distribution of individual amino acids
    '''Takes AA sequences and returns a (100,) np.array of distribution features for all AAs (present or not)'''
    L = len(X)
    result = []
    for aa in AAs:
        aas = [pos for pos,char in enumerate(a) if char==aa]
        l = len(aas)
        if l>2:
            t = np.ediff1d(aas)
            tm = t.mean()
            tsd = t.std()
            tp1 = np.percentile(t,5)
            tp2 = np.percentile(t,95)
            return np.hstack([l/L,tm,tsd,tp1,tp2])
        else:
            return np.hstack([l/L,0,0,0,0])


def statpack1(X):
    '''Takes numerical array and returns a (8,) np.array of distribution features'''
    m = X.mean()
    sd = X.std()
    p1 = np.percentile(X,5)
    p2 = np.percentile(X,95)
    t = np.ediff1d(X)
    tm = t.mean()
    tsd = t.std()
    tp1 = np.percentile(t,5)
    tp2 = np.percentile(t,95)
    return np.hstack([m,sd,p1,p2,tm,tsd,tp1,tp2])


def statpack2(X):
    '''Takes 2d numerical array and returns a (64,) np.array of distribution features'''
    arr = {}
    for z in [0,1]:
        m = statpack1(X.mean(axis=z))
        sd = statpack1(X.std(axis=z))
        p1 = statpack1(np.percentile(X,5,axis=z))
        p2 = statpack1(np.percentile(X,95,axis=z))
        arr[z] = np.hstack([m,sd,p1,p2])
    return np.hstack(arr.values())