package ModelFinal;

import HAL.GridsAndAgents.AgentGrid2D;
import HAL.GridsAndAgents.AgentSQ2Dunstackable;
import HAL.GridsAndAgents.PDEGrid2D;
import HAL.Gui.GridWindow;
import HAL.Rand;
import HAL.Tools.FileIO;
import HAL.Tools.Internal.Gaussian;
import HAL.Util;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.*;
import java.time.LocalDateTime; // import the LocalDateTime class

import static HAL.Util.*;

class ExampleCell extends AgentSQ2Dunstackable<ExampleGrid> {
    // declare individual cell attributes
    public int cellID; // cell ID
    public int type; // cell type
    public double IMTAge; // cell age
    public double cancerVariedIMT; // cancer cell inter mitotic age
    public double hd; // cancer cell threshold for death
    public double hp; // cancer cell threshold for proliferation
    public double activationProbability; // stroma cell probability of reactive phenotype
    public int activatedStromaNeighbours; // number of activated stroma neighbours of cancer cell

    // declare parameters

    // grid parameters
    static int x = 300; // number of grid locations horizontal

    static int y = 300; // number of grid locations vertical
    static double xDomain = 0.3; // (cm) horizontal length of sample space
    static double yDomain = 0.3; // (cm) vertical length of sample space
    static double deltaX = xDomain / x; // (cm) length of horizontal gridpoint
    static double deltaY = yDomain / y; // (cm) length of vertical gridpoint
    static double spaceConversion = (1 / deltaX) * (1 / deltaY); // space conversion parameter (used for PDEGrids)

    // time parameters
    static double timestepHour = 0.044; // (hours) length of timestep
    static double timestepPerDay = 24 / timestepHour; // number of timesteps in day
    static double timestepDay = timestepHour / 24;

    // agent paramenters
    static int STROMA = -1; // set code for stroma cells
    static int VESSEL = -2; // set code for vessel cells
    static int CANCER = -3; // set code for cancer cells
    static int ACTIVATEDSTROMA = -4; // set code for activated stroma cells
    static int QUIESCENTCANCER = -5; // set code for quiescent cancer cells
    static int UNREACTIVESTROMA = -6; // set code for unreactive stroma

    // vessel parameters
    static double sigmaMean = 0.016; // (cm) average distance between vessels
    static int nV = (int) ((xDomain * yDomain) / Math.pow(sigmaMean, 2)); // number of vessels in tissue

    // stroma parameters
    static int stromaNCI = 6; // stroma contact inhibition (number of empty neighbour cells required for division)

    static double stromaDieProbDays = 0.0042; // (days) stroma turnover probability
    static double stromaDieProb = stromaDieProbDays / timestepPerDay; // stroma turnover probability
    static int stromaIMTDays = 1; // (days) stroma cell inter mitotic time
    static int stromaIMT = stromaIMTDays * (int)timestepPerDay; // (timesteps) stroma cell inter mitotic time
    static double stromaDrugThresholdReactivity = 0.93; // stroma drug threshold reactivity
    static double stromaActivationNormalisationProbability = 10 * stromaDieProb;
    // stroma activation normalisation probability
    static double reactiveStromaProportion = 0.5; //proportion of stroma cell that can be activated


    // cancer parameters
    static int cellIDCount = 0;

    static int cancerIMTDays = 1; // (days) cancer cell inter mitotic time
    static int cancerIMT = cancerIMTDays *(int)timestepPerDay; // (timesteps) cancer cell inter mitotic time
    static boolean cancerEliminated = false; // boolean for cancer eliminated

    // drug parameters
    static double drugDiffusionCoefficientDays = 1E-5; // (cm^2 per day) diffusion coefficient drug
    static double drugDiffusionCoefficientTimestep = spaceConversion * drugDiffusionCoefficientDays * timestepDay;
    // (scaled) diffusion coefficient drug
    static double drugRemovalRateVesselDays = 500; // (days) removal rate of drug at vessel site
    static double drugRemovalRateVesselTimestep = drugRemovalRateVesselDays / timestepPerDay;
    // (timesteps) removal rate of drug at vessel site
    static double drugConcentrationVessel = 1; // drug concentration on delivery

    // proliferation signal parameters
    static double autocrineProliferationSignalProductionDays = 1.25; // (days) autocrine proliferation signal rate
    static double paracrineProliferationSignalProductionDays = 4.242; // (days) paracrine proliferation signal rate
    static double proliferationDegradationDrugDays = 6.68; // (days) proliferation signal degradation rate
    static double autocrineProliferationSignalProductionTimestep = autocrineProliferationSignalProductionDays
            / timestepPerDay; // (timesteps) autocrine proliferation signal rate
    static double paracrineProliferationSignalProductionTimestep = paracrineProliferationSignalProductionDays
            / timestepPerDay; // (timesteps) paracrine proliferation signal rate
    static double proliferationDegradationDrugTimestep = proliferationDegradationDrugDays / timestepPerDay;
    // (timesteps) proliferation signal degradation rate
    static double proliferationCancerSite = 0.256; // proliferation signal at initial cancer cell site


    // treatment parameters
    static int holidayDaysDrug = 20; // number of holiday days for drug


    public void stromaDieCell(double stromaDieProb) {
        // turnover of stroma cells with given probability
        if (G.rng.Double() < stromaDieProb) {
            //cell will die
            Dispose();
        }
    }

    public void stromaDivide(int position, int coordI) {
        // dividing cell dies and new cell of same type is created in same position
        int cellType = this.type;
        Dispose();
        ExampleCell resetCell = G.NewAgentSQ(coordI);
        resetCell.type = cellType;
        resetCell.IMTAge = 0;
        // create daughter cell and assign attributes
        ExampleCell daughterCell = G.NewAgentSQ(position);
        daughterCell.type = cellType;
        daughterCell.IMTAge = 0;
    }

    public void divideStromaCell(int stromaNCI) {
        // get x and y coordinates of cell
        int coordI = this.Isq();
        int options = MapEmptyHood(G.divHood);
        if (options >= stromaNCI) {
            int position = G.divHood[G.rng.Int(options)];
            stromaDivide(position, coordI);
        }
    }

    public boolean checkCancerNeighbour(int CANCER, int QUIESCENTCANCER) {
        boolean cancerPresent = false;
        int [] neighbourhood = G.divHood;
        int hoodsize = MapOccupiedHood(neighbourhood);
        for (int i = 0; i < hoodsize; i++) {
            ExampleCell checkCell = G.GetAgent(neighbourhood[i]);
            if (checkCell.type == CANCER || checkCell.type == QUIESCENTCANCER) {
                cancerPresent = true;
            }
        }
        return cancerPresent;
    }

    public int checkActivatedStromaNeighbour(int ACTIVATEDSTROMA) {
        int activatedStromaCn = 0;
        int [] neighbourhood = G.divHood;
        int hoodsize = MapOccupiedHood(neighbourhood);
        for (int i = 0; i < hoodsize; i++) {
            ExampleCell checkCell = G.GetAgent(neighbourhood[i]);
            if (checkCell.type == ACTIVATEDSTROMA) {
                activatedStromaCn++;
            }
        }
        return activatedStromaCn;
    }

    public void cancerDieCell() {
        //cell will die
        Dispose();
    }

    public int cancerDivide(int position, int coordI, double cancerIMT, double proliferationCancerSite,
                            int countCellID) {
        // dividing cell dies and new cell of same type is created in same position
        int cellType = this.type;
        Dispose();
        ExampleCell resetCell = G.NewAgentSQ(coordI);
        resetCell.cellID = cellIDCount;
        countCellID++;
        resetCell.type = cellType;
        resetCell.IMTAge = 0;

        // no mutation
        resetCell.hd = 0.2;
        resetCell.hp = 0.8;
        resetCell.cancerVariedIMT =  cancerIMT * (1 + (2*Math.random() - 1)*0.1);

        // create daughter cell and assign attributes
        ExampleCell daughterCell = G.NewAgentSQ(position);
        daughterCell.cellID = cellIDCount;
        countCellID++;
        daughterCell.type = QUIESCENTCANCER;
        daughterCell.IMTAge = 0;

        // no mutation
        daughterCell.hd = 0.2;
        daughterCell.hp = 0.8;
        daughterCell.cancerVariedIMT =  cancerIMT * (1 + (2*Math.random() - 1)*0.1);

        G.proliferation.Add(daughterCell.Xsq(), daughterCell.Ysq(), proliferationCancerSite);
        G.proliferation.Update();

        return countCellID;
    }

    public int divideCancerCell(double cancerIMT, double proliferationCancerSite, int cancerCellBirthCounter) {
        // declare neighbourhood integers
        int options;
        // get x and y coordinates of cell
        int coordI = this.Isq();
        options = MapEmptyHood(G.divHood);
        if (options > 0) {
            int position = G.divHood[G.rng.Int(options)];
            cancerCellBirthCounter = cancerDivide(position, coordI, cancerIMT, proliferationCancerSite,
                    cancerCellBirthCounter);
        }
        return cancerCellBirthCounter;
    }
}

public class ExampleGrid extends AgentGrid2D<ExampleCell> {
    // random number
    Rand rng = new Rand();
    // neighbourhood function
    int [] divHood = Util.MooreHood(false);

    // declare PDE grids
    PDEGrid2D drug;
    PDEGrid2D proliferation;

    public ExampleGrid(int x, int y) {
        // link agents and grid
        super(x, y, ExampleCell.class);
        // create PDE grid for drug
        drug = new PDEGrid2D(x, y);
        // create PDE grid for proliferation signal
        proliferation = new PDEGrid2D(x, y);
    }

    public int [] treatmentScheduleTimesteps(double timestepPerDay, int treatmentDaysDrug, int holidayDaysDrug) {
        int treatmentTimesteps = (int) (treatmentDaysDrug * timestepPerDay); // number of treatment timesteps for drug
        int holidayTimesteps = (int) (holidayDaysDrug* timestepPerDay); // number of holiday timesteps for drug
        int [] schedule = {treatmentTimesteps, holidayTimesteps};
        return schedule;
    }

    public void ageCells(int STROMA, int ACTIVATEDSTROMA, int CANCER, int UNREACTIVESTROMA) {
        ShuffleAgents(rng);
        // age stroma cell
        for (ExampleCell cell : this) {
            if (cell.type == ACTIVATEDSTROMA || cell.type == STROMA || cell.type == CANCER ||
                    cell.type == UNREACTIVESTROMA) {
                double currentAge = cell.IMTAge;
                cell.IMTAge = currentAge + 1;
            }
        }
    }

    public void stromaDieCells(double stromaDieProb, int STROMA, int ACTIVATEDSTROMA, int UNREACTIVESTROMA){
        ShuffleAgents(rng);
        // check for stroma cell type and turnover
        for (ExampleCell cell:this) {
            if ((cell.type == STROMA || cell.type == ACTIVATEDSTROMA || cell.type == UNREACTIVESTROMA)) {
                cell.stromaDieCell(stromaDieProb);
            }
        }
    }

    public void stromaStatusUpdate(int STROMA, int ACTIVATEDSTROMA, double stromaDrugThresholdReactivity,
                                   int CANCER, int QUIESCENTCANCER, double stromaActivationNormalisationProbability) {
        ShuffleAgents(rng);
        // check for stroma cell neighbouring cancer cells and drug threshold
        for (ExampleCell cell:this) {
            boolean cancerNeighbour = false;
            // activate stroma if it has cancer neighbour and drug is equal or above threshold
            cancerNeighbour = cell.checkCancerNeighbour(CANCER, QUIESCENTCANCER);
            if (cell.type == STROMA && cancerNeighbour == true) {
                // remove test for drug concentration so that it is alway met.
                double drugAmount = drug.Get(cell.Isq());
                if (drugAmount >= stromaDrugThresholdReactivity && rng.Double() <
                        stromaActivationNormalisationProbability) {
                    cell.type = ACTIVATEDSTROMA;
                }
                // deactivate stroma if drug is under threshold
            } else if (cell.type == ACTIVATEDSTROMA) {
                double drugAmount = drug.Get(cell.Isq());
                if (drugAmount < stromaDrugThresholdReactivity) {
                    cell.type = STROMA;
                }
            }
        }
    }

    public void initialiseDrugValues(int VESSEL, double drugConcentrationVessel) {
        // initialise drug from vessel sites
        for (ExampleCell cell:this) {
            if (cell.type == VESSEL) {
                drug.Set(cell.Xsq(), cell.Ysq(), drugConcentrationVessel);
            }
        }
        drug.Update();
    }

    public int heavysideValue(double inputValue) {
        int outputValue;
        if (inputValue > 0) {
            outputValue = 1;
        } else {
            outputValue = 0;
        }
        return outputValue;
    }

    public void updatePDEValues(int x, int y, int timeCounter, int treatmentTimestepDrug, int holidayTimestepDrug,
                                int VESSEL, double drugDiffusionCoefficientTimestep, double drugConcentrationVessel,
                                double drugRemovalRateVesselTimestep, int CANCER, int QUIESCENTCANCER,
                                double autocrineProliferationSignalProductionTimestep, int ACTIVATEDSTROMA,
                                double paracrineProliferationSignalProductionTimestep,
                                double proliferationDegradationDrugTimestep, int multiplierDeltaTHat) {

        boolean treatmentOn = false;
        // check if treatment is on
        if (timeCounter % (treatmentTimestepDrug + holidayTimestepDrug) < treatmentTimestepDrug) {
            treatmentOn = true;
        }

        // apply diffusion to drug according to treatment schedule
        if (treatmentOn == true) {
            initialiseDrugValues(VESSEL, drugConcentrationVessel);
        }
        drug.Diffusion(drugDiffusionCoefficientTimestep);

        // remove drug at vessel sites
        for (ExampleCell cell : this) {
            if (cell.type == VESSEL) {
                drug.Add(cell.Isq(), -drugRemovalRateVesselTimestep * drug.Get(cell.Isq()));
            }
        }

        if (timeCounter % multiplierDeltaTHat == 0) {
            // autocrine and paracrine increase
            for (ExampleCell cell : this) {
                if (cell.type == CANCER || cell.type == QUIESCENTCANCER) {
                    int valueHeavyside = heavysideValue(1 - proliferation.Get(cell.Isq()));
                    int activatedStromaNeighbourCn = cell.checkActivatedStromaNeighbour(ACTIVATEDSTROMA);
                    double paracrineInput = activatedStromaNeighbourCn * paracrineProliferationSignalProductionTimestep
                            * multiplierDeltaTHat;
                    double proliferationToAdd = (autocrineProliferationSignalProductionTimestep * multiplierDeltaTHat +
                            paracrineInput) * valueHeavyside;
                    proliferation.Add(cell.Isq(), proliferationToAdd);
                }
            }

            // degradation due to drug
            for (int i = 0; i < x * y; i++) {
                double drugDegradationValue = -proliferationDegradationDrugTimestep * multiplierDeltaTHat * drug.Get(i)
                        * proliferation.Get(i);
                proliferation.Add(i, drugDegradationValue);
            }
            // update proliferation grid
            proliferation.Update();
        }

        // update drug grid
        drug.Update();
    }

    public double drugMeanField() {
        double mFValue = drug.GetAvg();
        return mFValue;
    }

    public void cancerDieProliferationUpdate(int CANCER, int QUIESCENTCANCER) {
        ShuffleAgents(rng);
        for (ExampleCell cell:this) {
            if (cell.type == CANCER || cell.type == QUIESCENTCANCER) {
                double proliferationAmount = proliferation.Get(cell.Isq());
                if (proliferationAmount < cell.hd) {
                    // cancer cell dies
                    cell.cancerDieCell();
                    //System.out.println("Cancer Dies");
                } else if (proliferationAmount < cell.hp) {
                    // cancer cell is quiescent
                    cell.type = QUIESCENTCANCER;
                } else if (proliferationAmount >= cell.hp) {
                    // cancer cells is not quiescent
                    cell.type = CANCER;
                }
            }
        }
    }

    public int divideCells(int stromaNCI, int STROMA, int ACTIVATEDSTROMA, int stromaIMT, int CANCER, double cancerIMT,
                           double proliferationCancerSite, int UNREACTIVESTROMA, int cellIdentCount){
        ShuffleAgents(rng);
        // check cell type and its IMTAge and divide
        for (ExampleCell cell:this) {
            if ((cell.type == STROMA || cell.type == ACTIVATEDSTROMA || cell.type == UNREACTIVESTROMA) && cell.IMTAge
                    >= stromaIMT) {
                cell.divideStromaCell(stromaNCI);
            } else if (cell.type == CANCER && cell.IMTAge >= cell.cancerVariedIMT) {
                cellIdentCount = cell.divideCancerCell(cancerIMT, proliferationCancerSite, cellIdentCount);
            }
        }
        return cellIdentCount;
    }

    public boolean testElimination(int CANCER, int QUIESCENTCANCER) {
        int cancerCn = 0;
        boolean cancerEliminated = false;
        for (ExampleCell cell:this) {
            if (cell.type == CANCER || cell.type == QUIESCENTCANCER) {
                cancerCn++;
            }
        }
        if (cancerCn == 0) {
            System.out.println("cancer eliminated!");
            cancerEliminated = true;
        }
        return cancerEliminated;
    }
    public int countInitialStroma(String initialStromaCn) throws IOException {
        // Get initial stroma count from homeostatic run
        // Reads file and stores as string
        FileReader file = new FileReader(new File(initialStromaCn));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = br.readLine();
        // Tokenise each element in string
        StringTokenizer defaultTokenizer = new StringTokenizer(temp);
        // Store each token as element in double array
        int setUp;
        int convertedValue = Integer.parseInt(defaultTokenizer.nextToken());
        setUp = convertedValue;
        return setUp;
    }

    public double[][] initialStromaSetUp(String initialStromaPositions, int initialStromaCount) throws IOException {
        // Get initial stroma sites from homeostasis run
        // Reads file and stores as string
        FileReader file = new FileReader(new File(initialStromaPositions));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = "";
        double[][] setUp = new double[initialStromaCount][2];
        int i = 0;
        int j = 0;
        while ((temp = br.readLine()) != null) {
            // Tokenise each element in string
            StringTokenizer defaultTokenizer = new StringTokenizer(temp, ",");
            // Store each token as element in double array
            while (defaultTokenizer.hasMoreTokens()) {
                double convertedValue = Double.parseDouble(defaultTokenizer.nextToken());
                setUp[i][j] = convertedValue;
                j++;
            }
            j = 0;
            i++;
        }
        return setUp;
    }

    public void initialiseStroma(double[][] homeostaticStroma, int STROMA, int UNREACTIVESTROMA,
                                 double reactiveStromaProportion) {
        for (int i = 0; i < homeostaticStroma.length; i++) {
            // create initial stroma cell and assign attributes
            ExampleCell newCell = NewAgentSQ((int)homeostaticStroma[i][0]);
            newCell.activationProbability = Math.random();
            if (newCell.activationProbability <= reactiveStromaProportion) {
                newCell.type = STROMA;
            } else {
                newCell.type = UNREACTIVESTROMA;
            }
            newCell.IMTAge = homeostaticStroma[i][1];
        }
    }

    public int countInitialCancer(String initialCancerCn) throws IOException {
        // Get initial stroma count from homeostatic run
        // Reads file and stores as string
        FileReader file = new FileReader(new File(initialCancerCn));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = br.readLine();
        // Tokenise each element in string
        StringTokenizer defaultTokenizer = new StringTokenizer(temp);
        // Store each token as element in double array
        int setUp;
        int convertedValue = Integer.parseInt(defaultTokenizer.nextToken());
        setUp = convertedValue;
        return setUp;
    }

    public double[][] initialCancerSetUp(String initialCancerPositions, int initialCancerCount) throws IOException {
        // Get initial cancer sites from homeostasis run
        // Reads file and stores as string
        FileReader file = new FileReader(new File(initialCancerPositions));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = "";
        double[][] setUp = new double[initialCancerCount][6];
        int i = 0;
        int j = 0;
        while ((temp = br.readLine()) != null) {
            // Tokenise each element in string
            StringTokenizer defaultTokenizer = new StringTokenizer(temp, ",");
            // Store each token as element in double array
            while (defaultTokenizer.hasMoreTokens()) {
                double convertedValue = Double.parseDouble(defaultTokenizer.nextToken());
                setUp[i][j] = convertedValue;
                j++;
            }
            j = 0;
            i++;
        }
        return setUp;
    }

    public double[] initialProliferationSignal(String initialProliferationValues, int x, int y) throws IOException {
        // Get initial cancer sites from homeostasis run
        // Reads file and stores as string
        FileReader file = new FileReader(new File(initialProliferationValues));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = "";
        double[] setUp = new double[x * y];
        int i = 0;
        while ((temp = br.readLine()) != null) {
            // Tokenise each element in string
            StringTokenizer defaultTokenizer = new StringTokenizer(temp, "/n");
            // Store each token as element in double array
            while (defaultTokenizer.hasMoreTokens()) {
                double convertedValue = Double.parseDouble(defaultTokenizer.nextToken());
                setUp[i] = convertedValue;
            }
            i++;
        }
        return setUp;
    }

    public void meanFieldDrug(double[] timeVector,double[] meanFieldDrug, FileIO MFDrug) {
        MFDrug.Write("time,drug MF\n");
        for (int i = 0; i < meanFieldDrug.length; i++) {
            MFDrug.Write(timeVector[i] + "," + meanFieldDrug[i] + "\n");
        }
    }

    public void populationData(int simNo, double[] timeVector, double[] stromaPopulation, double[] cancerPopulation,
                               double[] activatedStromaPopulation, double[] quiescentCancerPopulation,
                               FileIO populations) {
        for (int i = 0; i < timeVector.length; i++) {
            populations.Write(simNo + "," + timeVector[i] + "," + stromaPopulation[i] + "," + cancerPopulation[i]
                    + "," + activatedStromaPopulation[i] + "," + quiescentCancerPopulation[i] + "\n");
        }
    }

    public int initialiseCancer(double[][] initialCancer, int CANCER, int QUIESCENTCANCER, int cellCountID) {
        for (int i = 0; i < initialCancer.length; i++) {
            // create initial stroma cell and assign attributes
            ExampleCell newCell = NewAgentSQ((int)initialCancer[i][0]);
            newCell.cellID = cellCountID;
            cellCountID++;
            if ((int)initialCancer[i][1] == 11) {
                newCell.type = CANCER;
            } else if ((int)initialCancer[i][1] == 22) {
                newCell.type = QUIESCENTCANCER;
            }
            newCell.IMTAge = initialCancer[i][2];
            newCell.cancerVariedIMT = initialCancer[i][3];
            newCell.hd = 0.2;
            newCell.hp = 0.8;
        }
        return cellCountID;
    }

    public void neighboursActivatedStroma(int CANCER, int ACTIVATEDSTROMA) {
        for (ExampleCell cell:this) {
            if (cell.type == CANCER) {
                cell.activatedStromaNeighbours = cell.checkActivatedStromaNeighbour(ACTIVATEDSTROMA);
            }
        }
    }

    public void initialiseProliferationSignal(double[] proliferationSignalInit) {
        // update proliferation signal PDEGrid
        for (int i = 0; i < proliferationSignalInit.length; i++) {
            proliferation.Add(i, proliferationSignalInit[i]);
        }
        proliferation.Update();
    }

    public int [][] readVesselInitialisation(String initialVesselPositions, int nV) throws IOException {
        FileReader file = new FileReader(new File(initialVesselPositions));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = "";
        int[][] setUp = new int[nV][2];
        int i = 0;
        int j = 0;
        while ((temp = br.readLine()) != null) {
            // Tokenise each element in string
            StringTokenizer defaultTokenizer = new StringTokenizer(temp, ",");
            // Store each token as element in double array
            while (defaultTokenizer.hasMoreTokens()) {
                int convertedValue = Integer.parseInt(defaultTokenizer.nextToken());
                setUp[i][j] = convertedValue;
                j++;
            }
            j = 0;
            i++;
        }
        return setUp;
    }

    public void initialiseVesselCells(int[][] vesselSetUp, int nV, int VESSEL) {
        for (int j = 0; j < nV; j++) {
            NewAgentSQ(vesselSetUp[j][0], vesselSetUp[j][1]).type = VESSEL;
        }
    }

    public double populationType(int TYPE) {
        double typeCount = 0;
        for (ExampleCell cell:this) {
            if (cell.type == TYPE) {
                typeCount++;
            }
        }
    return typeCount;
    }

    public void positionData(int i, int j, int k,  int x, int y, String resultsFolder, int multiplierDeltaTHat) {
        // create csv files
        int fileLabel = j / multiplierDeltaTHat;
        FileIO positionInformation = new FileIO(resultsFolder + "PositionData_treatment_" + k + "_sim_" + i +
                "_timestep_" + fileLabel + ".csv", "w");
        // Write data headings
        positionInformation.Write("cell position,cell type,drug value,proliferation value\n");
        for (int l = 0; l < x * y; l++) {
            if (GetAgent(l) == null) {
                positionInformation.Write(l + "," + "0" + "," + drug.Get(l) + "," + proliferation.Get(l) + "\n");
            } else {
                positionInformation.Write(l + "," + GetAgent(l).type + "," + drug.Get(l) + "," +
                        proliferation.Get(l) + "\n");
            }
        }
        // close timestep data files
        positionInformation.Close();
    }

    public int neighbourTypeCount(int position, int cellType) {
        // counts the number of agents of cellType in neighbourhood
        int typeCount = 0;
        int neighbourhood = MapHood(divHood, position);
        for (int i = 0; i < neighbourhood; i++) {
            if (GetAgent(divHood[i]) != null) {
                if (GetAgent(divHood[i]).type == cellType) {
                    typeCount++;
                }
            }
        }
        return typeCount;
    }

    public void recordNeighbours(int x, int y, int i, int j, int k, int STROMA, int ACTIVATEDSTROMA,
                                 int UNREACTIVESTROMA, int CANCER, int QUIESCENTCANCER, int VESSEL,
                                 String resultsFolder, int multiplierDeltaTHat) {
        // create csv files
        int fileLabel = j / multiplierDeltaTHat;
        FileIO neighbourInformation = new FileIO(resultsFolder + "NeighbourData_treatment_" + k + "_sim_" + i
                + "_timestep_" + fileLabel + ".csv", "w");
        // Write data headings
        neighbourInformation.Write("cell position,cancer neighbours,quiescent cancer neighbours," +
                "activated stroma neighbours,vessel neighbours,stroma neighbours,unreactive stroma neighbours," +
                "empty neighbours\n");
        for (int l = 0; l < x * y; l++) {
            int cancerNeighbourCount = neighbourTypeCount(l, CANCER);
            int quiescentCancerNeighbourCount = neighbourTypeCount(l, QUIESCENTCANCER);
            int activatedStromaNeighbourCount = neighbourTypeCount(l, ACTIVATEDSTROMA);
            int vesselNeighbourCount = neighbourTypeCount(l, VESSEL);
            int stromaNeighbourCount = neighbourTypeCount(l, STROMA);
            int unreactiveStromaNeighbourCount = neighbourTypeCount(l, UNREACTIVESTROMA);
            int emptyNeighbourCount = MapEmptyHood(divHood, l);
            neighbourInformation.Write(l + "," + cancerNeighbourCount + "," + quiescentCancerNeighbourCount + ","
                    + activatedStromaNeighbourCount + "," + vesselNeighbourCount + "," + stromaNeighbourCount + "," +
                    unreactiveStromaNeighbourCount + "," + emptyNeighbourCount + "\n");
        }
        neighbourInformation.Close();
    }

}
