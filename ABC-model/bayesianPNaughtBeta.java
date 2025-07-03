package ModelFinal;

import HAL.Tools.FileIO;

import java.io.IOException;

public class bayesianPNaughtBeta extends ExampleCell {
    public static void main(String[]args) throws IOException {


        int noPriors = 1000;

        // name destination folder
        String destinationFolder = "";

        // log start time of algorithm
        long startTime = System.nanoTime();

        int noDays = 29;
        int timesteps = (int) ((noDays + 1) * timestepPerDay);

        double [][] bayesianPriors = readPriors(destinationFolder + "BayesianPriorsPNaughtBeta.csv", noPriors);

        int noSims = 10;

        FileIO cancerTimeSeries = new FileIO(destinationFolder + "cancerTimeSeries_joint.csv", "w");
        // Write data headings
        cancerTimeSeries.Write("sim,time,cellCount,diameter,pNaught,beta\n");


        for (int i = 0; i < bayesianPriors.length; i++) {

            for (int k = 0; k < noSims; k++) {

                double betaVary = bayesianPriors[i][1] * timestepDay;

                boolean boundaryReached = false;
                // declare time and population storage arrays
                double[] timeVector = new double[timesteps + 1];
                double[] cancerPopulation = new double[timesteps + 1];
                double[] quiescentCancerPopulation = new double[timesteps + 1];

                // set up model
                ExampleGrid model = new ExampleGrid(x, y);

                // read vessel initialisation
                int[][] vesselSetUp = model.readVesselInitialisation(destinationFolder +
                        "vesselLocations.csv", nV);

                // initialise vessel cells
                model.initialiseVesselCells(vesselSetUp, nV, VESSEL);

                // determine position and age of initial cancer population
                int initialCancerCount = model.countInitialCancer(destinationFolder +
                        "cancerInitialCount.csv");
                double[][] initialCancer = model.initialCancerSetUp(destinationFolder +
                        "locationsCancer.csv", initialCancerCount);

                // determine intial proliferation signal over domain
                double[] proliferationSignalInit = model.initialProliferationSignal(
                        destinationFolder + "valuesProliferationSignal.csv", x, y);

                // initialise proliferation signal
                model.initialiseProliferationSignal(proliferationSignalInit);

                // initialise cancer
                model.initialiseCancer(initialCancer, CANCER, QUIESCENTCANCER, cellIDCount, parentCellHp, parentCellHd);

                double averageDiameter = model.calculateAverageDiameter(x, y, CANCER, QUIESCENTCANCER);

                model.collectTimeSeriesDataBayesianPNaughtBeta(x, averageDiameter, k, 0 * timestepDay,
                        initialCancerCount, bayesianPriors[i][0], bayesianPriors[i][1], cancerTimeSeries);

                // main loop
                for (int j = 0; j < timesteps + 1; j++) {

                    // collect population data
                    double timeCount = j * 1 / timestepPerDay;
                    double cancerCn = model.populationType(CANCER);
                    double quiescentCancerCn = model.populationType(QUIESCENTCANCER);
                    timeVector[j] = timeCount;
                    cancerPopulation[j] = cancerCn;
                    quiescentCancerPopulation[j] = quiescentCancerCn;

                    // age cells
                    model.ageCells(CANCER);

                    // divide cells
                    model.divideCells(CANCER, cancerIMT, bayesianPriors[i][0], cellIDCount);

                    // cancer cell proliferation update and age cancer cells
                    model.cancerDieProliferationUpdate(CANCER, QUIESCENTCANCER);

                    // update PDEGrid values
                    model.updatePDEValuesNoDrug(CANCER, QUIESCENTCANCER, betaVary);

                    // increment timestep
                    model.IncTick();

                    // check for cancer elimination to break loop
                    cancerEliminated = model.testElimination(CANCER, QUIESCENTCANCER);

                    // check for cancer elimination
                    if (cancerEliminated == true) {
                        break;
                    }

                    // calculate average diameter
                    averageDiameter = model.calculateAverageDiameter(x, y, CANCER, QUIESCENTCANCER);

                    boundaryReached = model.collectTimeSeriesDataBayesianPNaughtBeta(x, averageDiameter, k,
                            (j + 1) * timestepDay, cancerCn + quiescentCancerCn,
                            bayesianPriors[i][0], bayesianPriors[i][1], cancerTimeSeries);

                    if (boundaryReached == true) {
                        break;
                    }


                }
            }
        }

        //populations.Close();
        cancerTimeSeries.Close();


        // log end time of algorithm
        long endTime = System.nanoTime();
        // time elapsed
        long timeElapsed = endTime - startTime;
        long timeSeconds = timeElapsed / 1000000000;
        long timeMinutes = timeSeconds / 60;

        // indicate simulation has finished
        System.out.println("Simulation finished!");
        System.out.println("Time taken: " + timeSeconds + " seconds");
        System.out.println("Time taken: " + timeMinutes + " minutes");
    }
}


