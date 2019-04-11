#ifndef SIMULATIONDEF
#define SIMULATIONDEF

#include <iostream>         // Input and output
#include <iomanip>          // Input and output formatting
#include <fstream>          // File streams

#include <armadillo>        // Matrix library

#include <random>           // Random numbers
#include <math.h>           // Mathmatical constants
#include <algorithm>        // Mathmatical constants

#include <vector>           // C++ standard vector
#include <string.h>         // Strings

#include <cassert>          // Assertions

#include <sys/types.h>      // Packages for the directory
#include <sys/stat.h>       //    information,
#include <dirent.h>         //    handling and etc.
#include <ctime>            // Time functions


/* Class to contain simulation parameters, and which drives the simulation forward */
class Simulation {
 private:

    int     N_0;                    //          The number of cells in the simulation intially
    double  P_0;                    // [1/µm^3] The density of invading phages in the simulation initially

    int     N_max;                  //          The maximum number of cells in the simulation
    int     M_max;                  //          The maximum number of phages in the simulation

    double  dt;                     // [hour]   Default size of the time step
    double  dT_c;                   // [hour]   Size of the time step which updates the cells
    double  dT;                     // [hour]   Size of the time step which updates everything else

    double  dR;                     // [µm]     Depth of the infected layer (for spawning)
    int     nSamp;                  //          Number of samples to save per simulation hour
    int     res;                    //          Resolution of the nutrient grid
    double  dGrid;                  //          Spacing of grid (in units of critical radius R)
    double  dShell;                 // [µm]     Thickness of agent layer (for biomass)

    double  n_0;                    // [1/µm^3] Initial concentration of nutrient
    double  g_max;                  // [1/hour] Maximal growth rate for the cells
    double  R;                      // [µm]     The length scale for division (Typical volume (0.6 - 0.7) µm^3)
    double  k;                      // [N]      Parameter for repulsive potential
    double  K;                      // [1/µm^3] Michaels-Menten constant for Monod growth
    double  gamma;                  //          Probability to infect cell (also contains adsorption)
    double  alpha;                  //          Probability for phage to go lysogenic
    int     beta;                   //          Multiplication factor phage
    double  delta;                  // [1/hour] Amount of phage decay
    double  epsilon;                // [1/div]  Probability for offspring to turn resistant
    double  r;                      //          Constant used in the time-delay mechanism
    double  T_i;                    // [hour]   Time when the phage infections begins (-1 disables phage infection)
    double  L_box;                  // [µm]     Length of boundary condition box

    int     M_tot;                  //          The accumulated number of phages
    int     numB;                   //          Current tally of suceptible cells
    int     numL;                   //          Current tally of lysogenic cells
    int     numI;                   //          Current tally of infected cells
    int     numD;                   //          Current tally of dead cells
    int     numR;                   //          Current tally of resistant cells

    double  eta;                    //          Amount of division noise (0.5+eta) * V and (0.5-eta) * V
    double  nu;                     //          Amount of displacement noise of the new poles (width of gaussian)

    double  D_B;                    // [µm^2/hour] Diffusion constant for the cells
    double  D_P;                    // [µm^2/hour] Diffusion constant for the phage
    double  D_n;                    // [µm^2/hour] Diffusion constant for the nutrient

    double  h_agar;                 // [µm]    Total height of the soft agar layers
    double  h_cell;                 // [µm]    Distance of the cell colony to the hard agar

    int     Time;                   // Integer to keep track of time

    int     maxStep;                // Maximum time-step used in adaptive algorithm
    double  sigma;                  // Maximum number of standard deviations to move in diffusions term
    int     nSkips;                 // Counts the number of time-steps skipped with adaptive algorithm
    int     posibleStep;            // Largest possible timestep at current time

    double  p;                      // Helper variable to control phage spawning

    int     debug;                  // Set the debug level
    bool    debugBool;
    bool    exit;                   // Exit bool; terminates execution
    bool    firstRun;               // Bool to indicate if this run is the first

    bool    nutrientField;          // Boolean to toggle between nutrient field and just nutrients

    bool    exportAny;              //
    bool    exportCellData;         //
    bool    exportColonySize;       // Booleans to control the export output
    bool    exportPhageData;        //
    bool    exportNutrient;         //

    bool    singleInfectedCell;     // Boolean to enable phages via a single infected cell
    bool    planarPhageInvasion;    // Boolean to enable phages invading in a plane
    bool    uniformPhageInvasion;   // Boolean to enable phages invading uniformly from entire space
    bool    manyInfectedCells;      // Boolean to enable phages invading from an infected surface and uniformly from entire space
    int     forcedThickness;        // "Boolean" to force an infected surface of a given number

    int     invasionType;           // Number to remember the chose invasion type.

    bool    reflectingBoundary;     // Boolean to enable the simple box boundary condition with reflective sides
    bool    absorbingBoundary;      // Boolean to enable the simple box boundary condition with absorbing sides
    bool    experimentBoundary;     // Boolean to enable the advanced petri dish boundary conditions

    std::string phageType;          // Contains the type of phage chosen (sets parameters according to De Paepe)

    double rngSeed;                 // The seed for the random number generator
    std::mt19937 rng;               // Mersenne twister, random number generator
    std::uniform_real_distribution  <double> rand;
    std::normal_distribution        <double> randn;

    std::ofstream f_cells;          // Filestream to save configuration of cells
    std::ofstream f_colonySize;     // Filestream to save size of colony (Volume, N, N_Lys, N_Lyt)
    std::ofstream f_phages;         // Filestream to save configuration of phages
    std::ofstream f_n;              // Filestream to save concentration of nutrient / amount of nutrient
    std::ofstream f_log;            // Filestream to save log.txt

    std::string path;               // Sets the path to store in

    double benchMark[20];           // Benchmarking vector (holds run-times for each function)

    // Coordinates of cells in the simulation
    std::vector< std::vector<double> > cells;

    // Coordinates of central biomass in the simulation
    std::vector<double> biomass;

    // Coordinates of phage in the simulation
    std::vector< std::vector<double> > phages;

    // Helper quantities to reduce computations
    std::vector<double> center; // The geometric center of the colony
    double r_max;               // The largest distance from cell to geometric center
    int I_max;                  // The ID of the cell furthest away from center

    // Concentration of nutrients in the simulation
    double nutrient;
    arma::cube nutrient_grid;

    // Clever Neighbourhood search grid
    std::vector< std::vector< std::vector< std::vector<unsigned short int> > > > CN;

    // Density of cells which can grow
    arma::icube B;

    // Discrete laplace operator
    arma::sp_mat lap;


 public:
    // Constructers
    explicit Simulation(int N_0);                               // Direct constructer
    explicit Simulation(const Simulation& other);               // Copy constructor

    // Driver
    int      Run(double T);                                     // Controls the evaluation of the simulation

 private:
    void     Initialize();                                      // Initialize the simulation
    void     CellUpdate();                                      // Updates the cells
    void     PhageUpdate();                                     // Updates the phages

    // Simulation functions
    void     CellMovement(int I, double *B);                    // Get movement vector for cell I
    void     PhageMovement(int I, double *P);                   // Get movement vector for phage I
    int      PhageInfection(int I);                             // Check for infection oppotunity for phage I

    void     GrowCell(int I);                                   // Update radius for a single cell (and grow new cells)
    bool     GrowBiomass();                                     // Update radius for the biomass and absorb/spawn new cells
    bool     GrowInfection(int I);                              // Update stage of infection, and create new phages.

    double   GrowthRate(double x, double y, double z);          // Returns the growthrate at the location
    void     NutrientGridUpdate();                              // Update the nutrient grid from diff. eqns.

    void     ApplyBoundaryConditions(int m);                    // Applies the boundary conditions to the m'th phage
    void     SpawnPhages();                                     // Spanws phages according to the spawning rules

    void     CreateCentralBiomass();                            // Replaces the most central bacteria with abstract biomass

    double   PotentialGradient(double d, double r1, double r2); // Returns the gradient of the potential

    double   FirstPassageTimeDistribution();                    // Returns the probabilty for a phage to enter computational region

    void     ComputeColonyExtent();                             // Computes the center and extent of the colony

    void     UpdateNearestNeighbourGrid();                      // Updates the CN grid if cells have been removed
    std::vector<int> NearestNeighbours(int i,int j,int k);      // Returns all cells in the vicinity of gridpoint (i,j,k)

 public:
    // Settings
    void     TimeStep(double dT);                               // Set the size of the time-step
    void     TimeStepSkip(int skip);                            // Set the maximum "skip" size
    void     MaxStep(double sigma);                             // Set the maximum diffusion length (in sigma)
    void     BoundaryType(int type);                            // Set the type of boundary conditios (1 = simple, 2 = experiment)
    void     GridSize(double dGrid);                            // Set the size of the clever neighbour grid
    void     ShellSize(double dShell);                          // Set the size of the shell margin when using central biomass
    void     MaxGrowthRate(double g_max);                       // Sets the maximum growthrate
    void     MaxColonySize(int N_max);                          // Sets the maximum colonysize

    void     PhageInvasionStartTime(double T_i);                // Sets the time when the phages should start infecting
    void     PhageInvasionType(int type);                       // Set the type of phage invasion (1 = single, 2 = planar, 3 = uniform, 4 = shell)
    void     PhageInvasionType(int type, int thickness);
    void     PhageInitialDensity(double P_0);                   // Sets initial density of the phages (1/µm^2)

    void     PhageType(std::string& phageType);                 // Sets the phage parameters to De Paepe values for "phageType"
    void     PhageDiffusionConstant(double D_P);                // Sets the diffusion constant of the phages
    void     PhageInfectionRate(double r);                      // Sets rate of the infection increaasing in stage
    void     PhageDecayRate(double delta);                      // Set the decay rate of the phages
    void     PhageLysogenFrequency(double alpha);               // Set the lysogen frequency of the phage
    void     PhageBurstSize(int beta);                          // Set the size of the bursts
    void     PhageAdsorptionParameter(double gamma);            // Changes the adsorption parameter gamma

    void     CellMutationRate(double epsilon);                  // Changes the mutation/resistance rate of the cells
    void     CellInitialDensity(double B_0);                    // Sets initial density of the bacteria (1/µm^2)


 private:
    // Helping functions
    int      RandI(int n);                                      // Returns random integer less than n
    double   Rand(double n);                                    // Returns random double less than n
    double   RandN(double m, double s);                         // Returns random normal dist. number with mean m and variance s^2

 public:
    int      NumberOfUninfectedCells();                         // Returns the current number of uninfected cells
    int      NumberOfPhages();                                  // Returns the current number of phages
    double   Biomass();                                         // Returns the total volume of all the cells
    double   GetGamma();                                        // Returns the current gamma value
    double   GetColonyExtent();                                 // Returns the current r_max
    void     SetRngSeed(int n);                                 // Sets the seed of the random number generator

 private:
    void     deb(int n);                                        // Debug function (prints "debug" and int)
    void     WriteLog(double T);                                // Write a log.txt file
    void     SaveState();                                       // Store the endstate of the simulation

 public:
    void     LoadState(std::string& loadPath);                  // Load a previous state

    void     Quiet();                                           // Set debug level to 0

    void     ClearErrors();                                     // Sets exit to false;

    void     SetSamples(int nSamp);                             // Set the number of output samples

    // File outputs
    void     ExportCellData();                                  //
    void     ExportColonySize();                                //
    void     ExportPhageData();                                 // Sets booleans for export functions
    void     ExportNutrient();                                  //

 private:
    void     ExportData(double t);                              // Master function to export the data
    void     f_ExportCellData(double t);                        // Export the position and size of the cells
    void     f_ExportColonySize(double t);                      // Export the volume of colony and number of cells.
    void     f_ExportPhageData(double t);                       // Export the position and size of the phages
    void     f_ExportNutrient(double t);                        // Export the amount of nutrient

    // Data handling
    void        OpenFileStream(std::ofstream& stream,           // Open filstream if not allready opened
                            std::string& fileName);
    std::string GeneratePath();                                 // Generates a save path for datafiles

 public:
    void        SetFolderNumber(int number);                    // Sets the folder number (useful when running parralel code)
    void        SetPath(std::string& path);                     // Sets the folder path (useful when running parralel code)

    // Get properties
    std::string GetPath();                                      // Returns the save path
    int         GetTime();                                      // Returns the time
    int         GetDeltaT();                                    // Returns the time-step dT

    // Operators
    Simulation& operator=(const Simulation& rhs);               // Copy assignment

    // Clean up
    void DeleteFolder();                                        // Delete the data folder
 private:
    void DeleteFolderTree(const char* directory_name);          // Delete folders recursively

 public:
    ~Simulation();                                              // Destructor

};

#endif

