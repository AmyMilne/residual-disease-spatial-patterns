ABC Model

All python code is written in Python version 3.9.13.

The model uses the HAL Java package available at https://halloworld.org/.

Download the following files and store in a destination folder.

    ABC_priors_joint.py
    ABC_model_joint.py
    ABC_priors_delta.py
    ABC_model_delta.py
    ExperimentalData.csv
    vesselLocations.csv
    cancerInitialCount.csv
    locationsCancer.csv
    valuesProliferationSignal.csv

Download the following zipped files and unzip into the HAL-master folder. Specify the destination folder in the code.

    growthABC.zip
    decayABC.zip
    validationABC.zip

Pipeline.

    1. Run ABC_priors_joint.py. This creates a csv of the joint prior distribution for pNaught and beta.
    2. Run growthABC/bayesianPNaughtBeta.java to create the data for the ABC model.
    3. Run ABC_model_joint.py. This creates a csv of the posterior distribtion for pNaught and beta.
    4. Run ABC_priors_delta.py. This creates a csv of the prior distribution for delta.
    5. Run decayABC/bayesianDelta.java to create the data for the ABC model. Input calibrated values for pNaught and beta.
    6. Run ABC_model_delta.py. This creates a csv of the posterior distribtion for delta.
    7. Run validationABC/testRegrowth.java with calibrated values for pNaught, beta and delta for validation.
