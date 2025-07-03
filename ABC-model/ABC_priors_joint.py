# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 13:53:23 2025

@author: amymm
"""

import numpy as np
import pandas as pd

# specify destination folder
destinationFolder = r''

#%% Priors - joint distribution for pNaught and beta

# random number generator
rng = np.random.default_rng()

# sample size
sampleSize = 1000
noSims = 10

# constraint data
hp = 0.8
minIMT = 0.5

# parameter X - p0
# parameter Y - beta
# parameter bounds as vectors

minXY = [0, 0]
maxXY = [hp, 4]

# generate random pair for prior joint distribution
randXY = np.round(rng.uniform(minXY, maxXY, size = (sampleSize,2)), 10)
randomSample_raw_joint = pd.DataFrame(randXY, columns=['p0', 'beta'])

# check constraint
randomSample_raw_joint['inSpace'] = np.where((randomSample_raw_joint['beta'] >= minIMT * (hp - randomSample_raw_joint['p0'])), 1, 0)

# filter prior joint distribution
randomSample_joint = randomSample_raw_joint[randomSample_raw_joint['inSpace'] == 1].copy()
priorData_joint = randomSample_joint[['p0', 'beta']].reset_index(drop=True).copy()

# store prior joint distribution as csv
priorData_joint.to_csv(destinationFolder + 'BayesianPriorsPNaughtBeta.csv', index = False, header = False)

