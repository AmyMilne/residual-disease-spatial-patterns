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

#%% Priors - joint distribution for pNaught and beta

# read prior joint distribution
priorData_joint = pd.read_csv(destinationFolder + 'BayesianPriorsPNaughtBeta.csv', index = False, header = False)

#%% ABC model - the data produced by the modified vehicle control ABM for 10 simulations of each prior pair is required in the input folder location

noSims = 10

# bayesian threshold for pNaught and beta parameter estimation
bayesianThreshold_joint = 0.018

# ABC using least squares
def bayesianTest(inputVector, testVector):
    leastSquares = np.nansum(np.square(inputVector - testVector)) / (len(inputVector))
    return leastSquares

# pNaught and beta priors
tempDF_priors_joint = priorData_joint

# import model output data
tempDF_data_joint = pd.read_csv(destinationFolder + 'cancerTimeSeries_joint.csv')

# create data collection array
collectDataArray_joint = []

# test each set of priors    
for k in range(0, len(tempDF_priors_joint)):
    # create least squares data collection array
    leastSquaresData_joint = np.zeros(noSims)
        
    # filter ABM data output for considered prior pair
    sampleData_raw_joint = tempDF_data_joint[(tempDF_data_joint['pNaught'] == tempDF_priors_joint.iloc[k,0]) & (tempDF_data_joint['beta'] == tempDF_priors_joint.iloc[k,1])].copy()
        
    # test for each simulation
    for i in range(0, noSims):
        # filter ABM data output for considered simulation
        tempData_joint = sampleData_raw_joint[sampleData_raw_joint['sim'] == i].copy()
        # normalise diameter to initial diameter
        tempData_joint['normalisedDiameter'] = tempData_joint['diameter'] / tempData_joint['diameter'].iloc[0]
        # create time in weeks data
        tempData_joint['weeks'] = tempData_joint['time'] / 7
            
        # create data array for time data used for ABC model
        sampleData_valuesRaw_joint = np.zeros(len(timeVector_vC))
        for l in range(0, len(sampleData_valuesRaw_joint)):
            tempData_reduced_joint = tempData_joint[tempData_joint['weeks'] >= timeVector_vC[l]].reset_index(drop=True).copy()
            if (tempData_reduced_joint.shape[0] != 0):
                sampleData_valuesRaw_joint[l] = tempData_reduced_joint['normalisedDiameter'].iloc[0]
            else:
                sampleData_valuesRaw_joint[l] = 'NaN'
        
        # ABC model for sample data for considered prior pair simulation
        leastSquaresData_joint[i] = bayesianTest(sampleData_valuesRaw_joint, vC_mean)
    
    # determine proportion of simulations below bayesian threshold and store data
    belowThreshold_joint = np.where(leastSquaresData_joint <= bayesianThreshold_joint, 1, 0).sum() / noSims
    collectDataArray_joint = np.append(collectDataArray_joint,[tempDF_priors_joint.iloc[k,0], tempDF_priors_joint.iloc[k,1], belowThreshold_joint])
    
# reshape ABC model data and store in dataframe
collectDataReshaped_joint = collectDataArray_joint.reshape(tempDF_priors_joint.shape[0], tempDF_priors_joint.shape[1] + 1)
DFoutput_joint = pd.DataFrame(data=collectDataReshaped_joint, columns=['p0', 'beta', 'bayesValue'])

# filter ABC model posterior data
DFoutput_joint_hits = DFoutput_joint[DFoutput_joint['bayesValue'] == 1]

# store ABC model posterior data as csv
DFoutput_joint_hits.to_csv(destinationFolder + 'pNaughtBetaPosteriors.csv', index = False)

