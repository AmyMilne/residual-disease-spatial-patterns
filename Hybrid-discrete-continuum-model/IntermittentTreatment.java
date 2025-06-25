package ModelFinal;

import HAL.Gui.GifMaker;
import HAL.Gui.GridWindow;
import HAL.Tools.FileIO;

import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;

public class IntermittentTreatment extends ExampleCell {
    public static void main(String[]args) throws IOException {

        // name results folder
        String inputFolder = ""; // specify input folder

        // name results folder
        String resultsFolder = ""; // specify output folder

        int deltaTHatMultiplier = 10;

        int noDays = 590;
        int timesteps = ((int) ((noDays + 1) * timestepPerDay));
        int timestepsHat = timesteps / 10;
        int [] treatmentVary = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 590};
        int noSims = 1;

        for (int k = 0; k < treatmentVary.length; k++) {

            // declare storage files
            FileIO populations = new FileIO(resultsFolder + "Populations_treatment_" + treatmentVary[k] +
                    ".csv", "w");
            populations.Write("sim,time,stroma,cancer,activated stroma,quiescent cancer\n");
            FileIO drugMF = new FileIO(resultsFolder + "DrugMF_treatment_" + treatmentVary[k] + ".csv",
                    "w");

            for (int i = 0; i < noSims; i++) {

                // declare time and population storage arrays
                double[] timeVector = new double[timestepsHat + 1];
                double[] stromaPopulation = new double[timestepsHat + 1];
                double[] cancerPopulation = new double[timestepsHat + 1];
                double[] activatedStromaPopulation = new double[timestepsHat + 1];
                double[] quiescentCancerPopulation = new double[timestepsHat + 1];
                double [] meanFieldDrug = new double[timestepsHat + 1];

                ExampleGrid model = new ExampleGrid(x, y);

                // read vessel initialisation
                int[][] vesselSetUp = model.readVesselInitialisation(inputFolder +
                        "vesselLocations.csv", nV);

                // initialise vessel cells
                model.initialiseVesselCells(vesselSetUp, nV, VESSEL);

                // determine position and age of initial stroma population
                int initialStromaCount = model.countInitialStroma(inputFolder + "stromaInitialCount.csv");
                double[][] initialStroma = model.initialStromaSetUp(inputFolder +
                        "locationsStroma.csv", initialStromaCount);

                // determine position and age of initial cancer population
                int initialCancerCount = model.countInitialCancer(inputFolder + "cancerInitialCount.csv");
                double[][] initialCancer = model.initialCancerSetUp(inputFolder +
                        "locationsCancer.csv", initialCancerCount);

                // determine intial proliferation signal over domain
                double[] proliferationSignalInit = model.initialProliferationSignal(inputFolder +
                        "valuesProliferationSignal.csv", x, y);

                // initialise stroma cells
                model.initialiseStroma(initialStroma, STROMA, UNREACTIVESTROMA, reactiveStromaProportion);

                // initialise proliferation signal
                model.initialiseProliferationSignal(proliferationSignalInit);

                // initialise cancer
                cellIDCount = model.initialiseCancer(initialCancer, CANCER, QUIESCENTCANCER, cellIDCount);

                // main loop
                for (int j = 0; j < timesteps + 1; j++) {
                    // create time counter
                    int timeCounter = j;

                    if (j % deltaTHatMultiplier == 0) {
                        int indexValue = j / deltaTHatMultiplier;
                        // collect population data
                        double timeCount = j * 1 / timestepPerDay;
                        double stromaCn = model.populationType(STROMA);
                        double cancerCn = model.populationType(CANCER);
                        double activatedStromaCn = model.populationType(ACTIVATEDSTROMA);
                        double quiescentCancerCn = model.populationType(QUIESCENTCANCER);
                        timeVector[indexValue] = timeCount;
                        stromaPopulation[indexValue] = stromaCn;
                        cancerPopulation[indexValue] = cancerCn;
                        activatedStromaPopulation[indexValue] = activatedStromaCn;
                        quiescentCancerPopulation[indexValue] = quiescentCancerCn;
                        meanFieldDrug[indexValue] = model.drugMeanField();
                    }

                    // store data - can be turned on if you would like to collect data.
                    // Note for one simulation of 590 days is approximately 1.2GB of data.
                    /*if ((j % 1000 == 0)) {
                        // record each position on domain data
                        model.positionData(i, j, k, x, y, runDateTime, resultsFolder, deltaTHatMultiplier);
                        model.recordNeighbours(x, y, i, j, k, STROMA, ACTIVATEDSTROMA, UNREACTIVESTROMA, CANCER
                        , QUIESCENTCANCER, VESSEL, resultsFolder, runDateTime, deltaTHatMultiplier);
                    }*/

                    if (j % deltaTHatMultiplier == deltaTHatMultiplier - 1) {
                        // age cells
                        model.ageCells(STROMA, ACTIVATEDSTROMA, CANCER, UNREACTIVESTROMA);

                        // divide cells
                        cellIDCount = model.divideCells(stromaNCI, STROMA, ACTIVATEDSTROMA, stromaIMT, CANCER,
                                cancerIMT / deltaTHatMultiplier, proliferationCancerSite, UNREACTIVESTROMA,
                                cellIDCount);

                        // stroma cells dying
                        model.stromaDieCells(stromaDieProb * deltaTHatMultiplier, STROMA, ACTIVATEDSTROMA,
                                UNREACTIVESTROMA);

                        // cancer cell proliferation update and age cancer cells
                        model.cancerDieProliferationUpdate(CANCER, QUIESCENTCANCER);

                        // stroma cells activated or normalised
                        model.stromaStatusUpdate(STROMA, ACTIVATEDSTROMA, stromaDrugThresholdReactivity, CANCER,
                                QUIESCENTCANCER,
                                stromaActivationNormalisationProbability
                                        * deltaTHatMultiplier);

                        // update activated stroma neighbours
                        model.neighboursActivatedStroma(CANCER, ACTIVATEDSTROMA);

                    }

                    int[] treatmentTimestep = model.treatmentScheduleTimesteps(timestepPerDay, treatmentVary[k],
                            holidayDaysDrug);


                    // update PDEGrid values
                    model.updatePDEValues(x, y, timeCounter, treatmentTimestep[0], treatmentTimestep[1], VESSEL,
                            drugDiffusionCoefficientTimestep, drugConcentrationVessel, drugRemovalRateVesselTimestep,
                            CANCER, QUIESCENTCANCER, autocrineProliferationSignalProductionTimestep,
                            ACTIVATEDSTROMA, paracrineProliferationSignalProductionTimestep,
                            proliferationDegradationDrugTimestep, deltaTHatMultiplier);

                    // increment timestep
                    model.IncTick();

                    // check for cancer elimination to break loop
                    cancerEliminated = model.testElimination(CANCER, QUIESCENTCANCER);

                    // check for cancer elimination
                    if (cancerEliminated == true) {
                        break;
                    }

                }

                // record population data as csv
                model.populationData(i, timeVector, stromaPopulation, cancerPopulation, activatedStromaPopulation,
                        quiescentCancerPopulation, populations);
                if (i == 0) {
                    model.meanFieldDrug(timeVector, meanFieldDrug, drugMF);
                }
            }
            populations.Close();
            drugMF.Close();
        }
    }
}


