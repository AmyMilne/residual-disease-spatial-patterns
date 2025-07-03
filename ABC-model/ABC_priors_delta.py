# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 13:53:23 2025

@author: amymm
"""

import numpy as np
import pandas as pd

# specify destination folder
destinationFolder = r''

#%% Priors - distribution for delta

# random number generator
rng = np.random.default_rng()

# sample size
sampleSize = 1000
noSims = 10

# parameter Z - delta

minZ = [6]
maxZ = [7]

# generate random prior distribution
randZ = np.round(rng.uniform(minZ, maxZ, size = (sampleSize,1)), 10)
priorData = pd.DataFrame(randZ, columns=['delta'])

# store prior joint distribution as csv
priorData.to_csv(destinationFolder + 'BayesianPriorsDelta.csv', index = False, header = False)

