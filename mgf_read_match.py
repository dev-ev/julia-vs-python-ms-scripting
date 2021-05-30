#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd


# Write a function for mgf reading into an object

# In[2]:


def load_mgf(fname):
    '''Read the file into one huge list, without pre-defined format'''
    FIELDS = ('TITLE=', 'RTINSECONDS=', 'PEPMASS=', 'CHARGE=', 'SCANS=')
    
    def format_precursor(spectrum):
        #Cover for a case when there's no precursor intensity
        if ' ' in spectrum['PEPMASS']:
            spectrum['PEPMASS'] = [
                float(x) for x in spectrum['PEPMASS'].split(' ')
            ]
        else:
            spectrum['PEPMASS'] = [ float(spectrum['PEPMASS']),]

        #Check the polarity, which may or mmay not be given after the digits
        # "2+", "3-" etc
        polarityMultiplier = 1
        if spectrum['CHARGE'][-1] == '-':
            polarityMultiplier = -1

        if not spectrum['CHARGE'][-1].isnumeric():
            spectrum['CHARGE'] = polarityMultiplier * int( spectrum['CHARGE'][:-1] )
        else:
            spectrum['CHARGE'] = int( spectrum['CHARGE'] )
        
        return True
    
    def ms_data_to_df(spectrum):

        spectrum['ms_data'] = pd.DataFrame(
            spectrum['ms_data'],
            columns = ('m/z', 'Intensity')
        )

        return True
    
    spectraList = []
    with open(fname, 'r') as fh:
        state = False
        for line in fh:
            if line[0].isnumeric() and state == True:
                spectrum['ms_data'].append(
                    [ float(x) for x in line.rstrip().split(' ') ]
                    )
            elif 'BEGIN IONS' in line:
                spectrum = {'ms_data': []}
                state = True
            elif 'END IONS' in line:
                #Do not add the spectrum to the list if it doesn't contain fragment masses
                if len(spectrum['ms_data']) > 0:

                    ms_data_to_df(spectrum)
                    format_precursor(spectrum)
                    
                    spectraList.append(spectrum)
                state = False
            else:
                for fieldName in FIELDS:
                    if fieldName in line and state == True:
                        spectrum[fieldName[:-1] ] = line.rstrip().split(fieldName)[1]
                
    return spectraList


# In[3]:


fname = 'Yeast_1000spectra.mgf'


# In[4]:


get_ipython().run_cell_magic('timeit', '-r 20', 'load_mgf(fname)')


# In[5]:


res = load_mgf(fname)
len(res)


# In[6]:


res[500]


# In[7]:


res[500]['ms_data']


# Let's now check each spectrum for known mass differences

# Monoisotopic masses of amino acids:

# In[8]:


AA_DELTAS = {
    'G': 57.02147, 'A': 71.03712, 'S': 87.03203, 'P': 97.05277, 'V': 99.06842, 
    'T': 101.04768, 'Ccam': 160.03065, 'Cmes': 148.996912, 'I/L': 113.08407,
    'N': 114.04293, 'D': 115.02695, 'Q': 128.05858, 'K': 128.09497, 'E': 129.0426,
    'M': 131.04049, 'Mox': 147.0354, 'H': 137.05891, 'F': 147.06842, 'R': 156.10112,
    'Y': 163.06333, 'W': 186.07932
}


# In[9]:


#Flatten the values from the dictionary
singleResDeltas = np.array(
    list( AA_DELTAS.values() ), dtype = 'float64'
)
print(singleResDeltas.dtype)
#Add doubly-charged and triply-charged mass Deltas (simply divide by 2 and 3)
singleResDeltas = np.concatenate(
    (
        singleResDeltas,
        singleResDeltas / 2,
        singleResDeltas / 3
    )
)
print(singleResDeltas.shape)
singleResDeltas[:5]


# Now take the spectra one-by-one, find pairwise mass differences and match them to the list.<br>
# * Calculate pairwise absolute differences between the 
# * Subtract the experimental mass Deltas from the theoretical
# * Calculate relative difference
# * Select the cases with the relative difference lower than threshold (matches)
# * Summarize and report the matches

# In[10]:


def find_matches(spectra, masses_to_match, rel_tolerance = 1e-5, float_arr_type = 'float64'):
    resDict = {
        'Spectrum_idx': np.array([], dtype='uint32'),
        'Exp_idx': np.array([], dtype='uint32'),
        'Library_idx': np.array([], dtype='uint32'),
        'Rel_error': np.array([], dtype=float_arr_type)
    }
    #Calculate the minimal value in the list for matching
    #and offset it by the matching tolerance
    minTheoVal = masses_to_match.min() * (1 - rel_tolerance)

    for idx, s in enumerate(spectra):
        #Calculate pairwise differences betweeen experimental values
        expDeltas = np.subtract.outer(
            s['ms_data']['m/z'].to_numpy(), s['ms_data']['m/z'].to_numpy()
        )            
        # Disregard relative deltas that are smaller than the lowest theoretical value
        expDeltas = expDeltas[ expDeltas > minTheoVal]
        # Calculate relative differences between experimental and theoretical values 
        relDeltasArr = np.divide(
            #Absolute values of the differences between masses
            np.abs(
                np.subtract.outer(
                    masses_to_match, expDeltas
                )
            ),
            #Means between the masses
            (np.add.outer(masses_to_match, expDeltas) / 2)
        )
        matchingInds = np.where(
            pd.DataFrame(
                relDeltasArr
            ).le(rel_tolerance) == True
        )
        numMatches = matchingInds[0].shape[0]
        if numMatches > 0:
            resDict['Spectrum_idx'] = np.append(
                resDict['Spectrum_idx'],
                np.array( [idx, ] * numMatches, dtype='uint32' )
            )
            resDict['Library_idx'] = np.append(
                resDict['Library_idx'], matchingInds[0]
            )
            resDict['Exp_idx'] = np.append(
                resDict['Exp_idx'], matchingInds[1]
            )
            resDict['Rel_error'] = np.append(
                resDict['Rel_error'],
                relDeltasArr[ matchingInds[0], matchingInds[1] ]
            )

    resDF = pd.DataFrame(resDict)
    return resDF


# In[11]:


get_ipython().run_cell_magic('timeit', '-r 5', 'find_matches(res, singleResDeltas, 1e-5)')


# In[12]:


matches = find_matches(res, singleResDeltas, rel_tolerance = 1e-5)
matches


# In[13]:


matches[ matches['Spectrum_idx'] == 1 ]


# In[ ]:




