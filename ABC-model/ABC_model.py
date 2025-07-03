# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 13:53:23 2025

@author: amymm
"""

import numpy as np
import pandas as pd

# specify input and output folders
inputFolder = r''
outputFolder = r''

#%% ExperimentalData - Experimental data excel spread sheet is required in the input folder location

# file path to experimental data
inputFile = inputFolder + r'ExperimentalData.xlsx'

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
priorData_joint.to_csv(outputFolder + 'BayesianPriorsPNaughtBeta.csv', index = False, header = False)

#%% ABC model - the data produced by the modified vehicle control ABM for 10 simulations of each prior pair is required in the input folder location

# bayesian threshold for pNaught and beta parameter estimation
bayesianThreshold_joint = 0.018

# ABC using least squares
def bayesianTest(inputVector, testVector):
    leastSquares = np.nansum(np.square(inputVector - testVector)) / (len(inputVector))
    return leastSquares

# pNaught and beta priors
tempDF_priors_joint = priorData_joint

# import model output data
tempDF_data_joint = pd.read_csv(inputFolder + 'cancerTimeSeries_joint.csv')

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
DFoutput_joint_hits.to_csv(outputFolder + 'pNaughtBetaPosteriors.csv', index = False)

#%% Priors - distribution for delta

# parameter Z - delta

minZ = [6]
maxZ = [7]

# generate random prior distribution
randZ = np.round(rng.uniform(minZ, maxZ, size = (sampleSize,1)), 10)
priorData = pd.DataFrame(randZ, columns=['delta'])

# store prior joint distribution as csv
priorData.to_csv(outputFolder + 'BayesianPriorsDelta.csv', index = False, header = False)

#%% ABC model - the data produced by the modified treatment ABM for 10 simulations of each prior is required in the input folder location

# bayesian threshold for delta parameter estimation
bayesianThreshold = 0.0005

# ABC using least squares
def bayesianTest(inputVector, testVector):
    leastSquares = np.nansum(np.square(inputVector - testVector)) / (len(inputVector))
    return leastSquares

# delta priors
tempDF_priors = priorData

# import model output data
tempDF_data = pd.read_csv(inputFolder + 'cancerTimeSeries.csv')

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
DFoutput_hits.to_csv(outputFolder + 'deltaPosteriors.csv', index = False)




