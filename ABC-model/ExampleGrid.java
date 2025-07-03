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
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;
import java.time.LocalDateTime; // import the LocalDateTime class

import static HAL.Util.*;

class ExampleCell extends AgentSQ2Dunstackable<ExampleGrid> {
    // declare individual cell attributes
    public int cellID;
    public int type;
    public double IMTAge;
    public double cancerVariedIMT;
    public double hd;
    public double hp;

    // declare parameters

    // grid parameters
    static int x = 500; // number of cells horizontal

    static int y = 500; // number of cells vertical
    static double xDomain = 0.5; // (cm) horizontal length of sample space
    static double yDomain = 0.5; // (cm) vertical length of sample space
    static double deltaX = xDomain / x; // (cm) length of horizontal gridpoint
    static double deltaY = yDomain / y; // (cm) length of vertical gridpoint
    static double spaceConversion = (1 / deltaX) * (1 / deltaY); // space conversion parameter (used for PDEGrids)

    // time parameters
    static double timestepHour = 0.44; // (hours) length of timestep
    static double timestepPerDay = 24 / timestepHour; // number of timesteps in day
    static double timestepDay = timestepHour / 24;

    // agent paramenters
    static int VESSEL = -2; // set code for vessel cells
    static int CANCER = -3; // set code for cancer cells
    static int QUIESCENTCANCER = -5; // set code for quiescent cancer cells

    // vessel parameters
    static double sigmaMean = 0.016; // (cm) average distance between vessels
    static int nV = (int) ((xDomain * yDomain) / Math.pow(sigmaMean, 2)); // number of vessels in tissue

    // cancer parameters
    static int cellIDCount = 0;

    // cancer proliferation signal threshold proliferation
    static double cancerIMTDays = 1; // (days) cancer cell inter mitotic time
    static double cancerIMT = cancerIMTDays *(int)timestepPerDay; // (timesteps) cancer cell inter mitotic time
    static boolean cancerEliminated = false; // boolean for cancer eliminated

    static double parentCellHd = 0.2;//
    static double parentCellHp = 0.8;//

    public void cancerDieCell() {
        //cell will die
        Dispose();
    }

    public int cancerDivide(int position, int coordI, double cancerIMT, double proliferationCancerSite,
                            int countCellID) {
        // dividing cell dies and new cell of same type is created in same position
        int cellType = this.type;
        double parentHd = this.hd;
        double parentHp = this.hp;
        Dispose();
        ExampleCell resetCell = G.NewAgentSQ(coordI);
        resetCell.cellID = countCellID;
        countCellID++;
        resetCell.type = cellType;
        resetCell.IMTAge = 0;


        // no mutation
        resetCell.hd = parentHd;
        resetCell.hp = parentHp;
        resetCell.cancerVariedIMT =  cancerIMT * (1 + (2*Math.random() - 1)*0.1);

        // create daughter cell and assign attributes
        ExampleCell daughterCell = G.NewAgentSQ(position);
        daughterCell.cellID = countCellID;
        countCellID++;
        daughterCell.type = QUIESCENTCANCER;
        daughterCell.IMTAge = 0;

        // no mutation
        daughterCell.hd = parentHd;
        daughterCell.hp = parentHp;
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
        options = MapEmptyHood(G.divHoodIncreased);
        if (options > 0) {
            int tempIndex = G.rng.Int(options);
            // increased divisionhood
            int position = G.divHoodIncreased[tempIndex];
            cancerCellBirthCounter = cancerDivide(position, coordI, cancerIMT, proliferationCancerSite,
                    cancerCellBirthCounter);
        }
        return cancerCellBirthCounter;
    }

    public static double [][] readPriors(String priorFile, int priorCount) throws IOException {
        double [][] tempPriors = new double [priorCount][2];
        FileReader file = new FileReader(new File(priorFile));
        // declare reading variables
        BufferedReader br = new BufferedReader(file);
        String temp = "";
        int i = 0;
        int j = 0;
        while ((temp = br.readLine()) != null) {
            // Tokenise each element in string
            StringTokenizer defaultTokenizer = new StringTokenizer(temp, ",");
            // Store each token as element in double array
            while (defaultTokenizer.hasMoreTokens()) {
                double convertedValue = Double.parseDouble(defaultTokenizer.nextToken());
                tempPriors[i][j] = convertedValue;
                j++;
            }
            j = 0;
            i++;
        }
        return tempPriors;
    }

}

public class ExampleGrid extends AgentGrid2D<ExampleCell> {
    // random number
    Rand rng = new Rand();
    // neighbourhood function
    int increasedDivisionhoodRadius = 6;
    int [] divHoodIncreased = Util.RectangleHood(false, increasedDivisionhoodRadius,
            increasedDivisionhoodRadius);

    // declare PDE grids
    PDEGrid2D proliferation;

    public ExampleGrid(int x, int y) {
        // link agents and grid
        super(x, y, ExampleCell.class);
        // create PDE grid for proliferation signal
        proliferation = new PDEGrid2D(x, y);
    }

    public void ageCells(int CANCER) {
        ShuffleAgents(rng);
        // age stroma cell
        for (ExampleCell cell : this) {
            if (cell.type == CANCER) {
                double currentAge = cell.IMTAge;
                cell.IMTAge = currentAge + 1;
            }
        }
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

    public void updatePDEValuesNoDrug(int CANCER, int QUIESCENTCANCER,
                                      double autocrineProliferationSignalProductionTimestep) {


        for (ExampleCell cell : this) {
            if (cell.type == CANCER || cell.type == QUIESCENTCANCER) {
                int valueHeavyside = heavysideValue(1 - proliferation.Get(cell.Isq()));
                double proliferationToAdd = autocrineProliferationSignalProductionTimestep * valueHeavyside;
                proliferation.Add(cell.Isq(), proliferationToAdd);
            }
        }
        // update proliferation grid
        proliferation.Update();
    }

    public int cancerDieProliferationUpdate(int CANCER, int QUIESCENTCANCER) {
        ShuffleAgents(rng);
        int cancerDieTemp = 0;
        for (ExampleCell cell:this) {
            if (cell.type == CANCER || cell.type == QUIESCENTCANCER) {
                double proliferationAmount = proliferation.Get(cell.Isq());
                if (proliferationAmount < cell.hd) {
                    // cancer cell dies
                    cell.cancerDieCell();
                    cancerDieTemp++;
                } else if (proliferationAmount < cell.hp) {
                    // cancer cell is quiescent
                    cell.type = QUIESCENTCANCER;
                } else if (proliferationAmount >= cell.hp) {
                    // cancer cells is not quiescent
                    cell.type = CANCER;
                }
            }
        }
        return cancerDieTemp;
    }

    public int divideCells(int CANCER, double cancerIMT, double proliferationCancerSite, int cellIdentCount){
        ShuffleAgents(rng);
        // check cell type and its IMTAge and divide
        for (ExampleCell cell:this) {
            if (cell.type == CANCER && cell.IMTAge >= cell.cancerVariedIMT) {
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
            cancerEliminated = true;
        }
        return cancerEliminated;
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

    public int initialiseCancer(double[][] initialCancer, int CANCER, int QUIESCENTCANCER, int cellCountID,
                                double parentCellHp, double parentCellHd) {
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
            newCell.hd = parentCellHd;
            newCell.hp = parentCellHp;
        }
        return cellCountID;
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

    public double calculateAverageDiameter(int x, int y, int CANCER, int QUIESCENTCANCER) {
        int tempMinX = x;
        int tempMaxX = 0;
        int tempMinY = y;
        int tempMaxY = 0;
        for (ExampleCell cell : this) {
            if (cell.type == CANCER || cell.type == QUIESCENTCANCER) {
                if (cell.Xsq() < tempMinX) {
                    tempMinX = cell.Xsq();
                }
                if (cell.Xsq() > tempMaxX) {
                    tempMaxX = cell.Xsq();
                }
                if (cell.Ysq() < tempMinY) {
                    tempMinY = cell.Ysq();
                }
                if (cell.Ysq() > tempMaxY) {
                    tempMaxY = cell.Ysq();
                }
            }
        }

        int xDiff = tempMaxX - tempMinX + 1;
        int yDiff = tempMaxY - tempMinY + 1;

        // note that this assumes a square domain
        double aveDiameter = 0.5 * (xDiff + yDiff);

        return aveDiameter;
    }

    public boolean collectTimeSeriesDataBayesianPNaughtBeta(int x, double aveDiameter, int simNoCurrent,
                                                            double currentTime, double countCancer,
                                                            double currentPNaught, double currentBeta,
                                                            FileIO resultsFile) {

        boolean tempBoundaryReached;
        if (aveDiameter != x) {
            tempBoundaryReached = false;
        } else {
            tempBoundaryReached = true;
        }

        resultsFile.Write(simNoCurrent + "," + currentTime + "," + countCancer + "," + aveDiameter + "," +
                currentPNaught + "," + currentBeta + "\n");
        return tempBoundaryReached;
    }

}
