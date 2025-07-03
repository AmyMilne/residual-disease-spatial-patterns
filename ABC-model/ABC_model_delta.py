# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 13:53:23 2025

@author: amymm
"""

import numpy as np
import pandas as pd

# specify destination folder
destinationFolder = r''

#%% ExperimentalData - Experimental data excel spread sheet is required in the input folder location

# file path to experimental data
inputFile = destinationFolder + r'ExperimentalData.xlsx'

# import experimental data
DFinputData_raw = pd.read_excel(inputFile, header = None)#, index_col=None, sheet_name = 'XenograftModelData')

# filter vehicle control and treatment data
vC = DFinputData_raw.iloc[3:7,1:11].to_numpy()
continuousTreatment = DFinputData_raw.iloc[3:, 11:].to_numpy()

# create data arrays for normalised experimental data
vC_norm = np.zeros(vC.shape)
continuousTreatment_norm = np.zeros(continuousTreatment.shape)

for i in range(0, vC.shape[1]):
    vC_norm[:,i] = vC[:,i] / vC[0,i]

for i in range(0, continuousTreatment.shape[1]):
    continuousTreatment_norm[:,i] = continuousTreatment[:,i] / continuousTreatment[0,i]

# mean diameter data across replicates
vC_mean = np.nanmean(vC_norm, axis = 1)
continuousTreatment_mean = np.nanmean(continuousTreatment_norm, axis = 1)[:3]

# create time vectors for growth and decay of diameters
timeVector_vC = np.array([0, 1, 2, 3])
timeVector_continuous = np.array([0, 1, 2])

#%% Priors - distribution for delta

# read prior joint distribution
priorData = pd.read_csv(destinationFolder + 'BayesianPriorsDelta.csv', index = False, header = False)

#%% ABC model - the data produced by the modified treatment ABM for 10 simulations of each prior is required in the input folder location

# bayesian threshold for delta parameter estimation
bayesianThreshold = 0.0005

noSims = 10

# ABC using least squares
def bayesianTest(inputVector, testVector):
    leastSquares = np.nansum(np.square(inputVector - testVector)) / (len(inputVector))
    return leastSquares

# delta priors
tempDF_priors = priorData

# import model output data
tempDF_data = pd.read_csv(destinationFolder + 'cancerTimeSeries.csv')

# create data collection array
collectDataArray = []

# test each prior    
for k in range(0, len(tempDF_priors)):
    # create least squares data collection array
    leastSquaresData = np.zeros(noSims)
        
    # filter ABM data output for considered prior
    sampleData_raw = tempDF_data[tempDF_data['delta'] == tempDF_priors.iloc[k]].copy()
        
    # test for each simulation
    for i in range(0, noSims):
        # filter ABM data output for considered simulation
        tempData = sampleData_raw[sampleData_raw['sim'] == i].copy()
        # normalise diameter to initial diameter
        tempData['normalisedDiameter'] = tempData['diameter'] / tempData['diameter'].iloc[0]
        # create time in weeks data
        tempData['weeks'] = tempData['time'] / 7
            
        # create data array for time data used for ABC model
        sampleData_valuesRaw = np.zeros(len(timeVector_continuous))
        for l in range(0, len(sampleData_valuesRaw)):
            tempData_reduced = tempData[tempData['weeks'] >= timeVector_continuous[l]].reset_index(drop=True).copy()
            if (tempData_reduced.shape[0] != 0):
                sampleData_valuesRaw[l] = tempData_reduced['normalisedDiameter'].iloc[0]
            else:
                sampleData_valuesRaw[l] = 'NaN'
        
        # ABC model for sample data for considered prior simulation
        leastSquaresData[i] = bayesianTest(sampleData_valuesRaw, continuousTreatment_mean)
    
    # determine proportion of simulations below bayesian threshold and store data
    belowThreshold = np.where(leastSquaresData <= bayesianThreshold, 1, 0).sum() / noSims
    collectDataArray = np.append(collectDataArray,[tempDF_priors.iloc[k,0], tempDF_priors.iloc[k,1], belowThreshold])
    
# reshape ABC model data and store in dataframe
collectDataReshaped = collectDataArray.reshape(tempDF_priors.shape[0], tempDF_priors.shape[1] + 1)
DFoutput = pd.DataFrame(data=collectDataReshaped, columns=['delta', 'bayesValue'])

# filter ABC model posterior data
DFoutput_hits = DFoutput[DFoutput['bayesValue'] == 1]

# store ABC model posterior data as csv
DFoutput_hits.to_csv(destinationFolder + 'deltaPosteriors.csv', index = False)




