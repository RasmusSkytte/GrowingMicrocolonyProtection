#include "Simulation.hpp"
#include <omp.h>

using namespace std;
using namespace arma;


// Constructers /////////////////////////////////////////////////////////////////////////
// Direct constructer
Simulation::Simulation(int N_0) {

    // Store the number of cells in the simulation
    this->N_0 = N_0;

    // Set some default parameters (initlize some default objects)
    P_0                 = 0;        // [1/µm^2] The density of invading phages in the simulation initially
    N_max               = 65536;    // Maximum number of cells in simulation (Unsigned short int)
    M_max               = 1e6;      // Maximum number of phages in simulation
    M_tot               = 0;        // Accumulated number of phages
    dt                  = 1e-6;     // [hour]   Default time step size
    dT_c                = dt;       // [hour]   Size of the time step which updates the phages
    dT                  = dt;       // [hour]   Size of the time step which updates everything else
    dR                  = 3.5;      // Depth of the infected layer (for spawning)
    nSamp               = 10;       // Number of samples to save per simulation hour
    dGrid               = 2;        // Spacing of grid (in units of critical radius R)
    dShell              = HUGE_VAL; // [µm]     Thickness of agent layer (for biomass)

    n_0                 = 0.1;      // [1/µm^3] Initial concentration of nutrient
    g_max               = 10.0/6.0; // [1/hour] Maximal growth rate for the cells
    R                   = 0.782;    // [µm]     The length scale for division (Typical volume 1.3 µm^3)
    k                   = 1e3;      // [N*m]    Parameter for repulsive potential
    K                   = 0.0;      //          Michaels-Menten FACTOR for Monod growth (K = 0.0 disables Monod behavior)
    gamma               = 1/dT;     //          Probability to infect cell
    alpha               = 0.00;     //          Probability for phage to go lysogenic
    beta                = 400;      //          Multiplication factor phage
    delta               = 0.003;    // [1/hour] Rate of phage decay
    epsilon             = 0;        //          Probability for offspring to turn resistant
    r                   = 10.0*2/3; //          Constant used in the time-delay mechanism
    T_i                 = -1;       // [hours]  Time when the phage infections begins (less than 0 disables phage infection)
    L_box               = 50;       // [µm]     Length of boundary condition box

    numB                = 0;        // Current tally of suceptible cells
    numL                = 0;        // Current tally of lysogenic cells
    numI                = 0;        // Current tally of infected cells
    numD                = 0;        // Current tally of dead cells
    numR                = 0;        // Current tally of resistant cells

    eta                 = 0.1;      // Amount of division noise (width of gaussian)
    nu                  = R/4;      // Amount of displacement noise during division. x = x0 + rand(nu) etc.

    D_B                 = 0.0;      // [µm^2/hour] Diffusion constant for the cells
    D_P                 = 13000;    // [µm^2/hour] Diffusion constant for the phage
    D_n                 = 4e5;      // [µm^2/hour] Diffusion constant for the nutrient

    h_agar              = 250;      // [µm]    Total height of the soft agar layers
    h_cell              = 125;      // [µm]    Distance of the cell colony to the hard agar

    p                   = 0;        // Helper variable to control phage spawning

    Time                = 0;        // Counter for how many timesteps have passed
    maxStep             = 1;        // Maximum time-step used in adaptive algorithm
    sigma               = HUGE_VAL; // Maximum number of standard deviations to move in diffusions term
    nSkips              = 0;        // Counts the number of time-steps skipped with adaptive algorithm
    posibleStep         = 1;        // Largest possible "step" at current time

    debug               = 1;        // The amount of information to print to terminal
    debugBool           = false;
    exit                = false;    // Boolean to control early exit

    firstRun            = true;     // Bool to indicate if this run is the first

    nutrientField       = false;    // Boolean to toggle between nutrient field and just nutrients

    exportAny           = false;    //
    exportCellData      = false;    //
    exportColonySize    = false;    // Booleans to control the export output
    exportPhageData     = false;    //
    exportNutrient      = false;    //

    singleInfectedCell  = false;    // Boolean to enable phages via a single infected cell
    planarPhageInvasion = false;    // Boolean to enable phages invading in a plane
    uniformPhageInvasion= false;    // Boolean to enable phages invading uniformly from entire space
    manyInfectedCells   = false;    // Boolean to enable surface infection and uniformly from entire space
    forcedThickness     = 0;        // "Boolean" to force an infected surface of a given number

    invasionType        = 0;        // Number to remember the invasion type

    reflectingBoundary  = false;    // Boolean to enable the simple box boundary condition with reflective sides
    absorbingBoundary   = true;     // Boolean to enable the simple box boundary condition with absorbing sides
    experimentBoundary  = false;    // Boolean to enable the advanced petri dish boundary conditions

    phageType = "default";          // Contains the type of phage chosen (sets parameters according to De Paepe)
    rngSeed = -1;                   // Random number seed  ( set to -1 if unused )

    // Initilize the benchmarking vector
    for (int j = 0; j < 20; j ++) {
        benchMark[j] = 0.0;
    }

    // Set the resolution of the grid to match Lbox
    res = ceil(L_box/(2*dGrid*R)-1/2);

    // Set the number of threads to use
    omp_set_num_threads(1);
};


// Copy constructor
Simulation::Simulation(const Simulation& other) {

    N_0                 = other.N_0;                        // Store the initial number of cells in the simulation
    P_0                 = other.P_0;                        // [1/µm^2] The density of invading phages in the simulation initially
    N_max               = other.N_max;                      // Maximum number of cells in simulation
    M_max               = other.M_max;                      // Maximum number of phages in simulation
    M_tot               = other.M_tot;                      // Accumulated number of phages
    dt                  = other.dt;                         // [hour]   Default time step size
    dT_c                = other.dT_c;                       // [hour]   Size of the time step which updates the cells
    dT                  = other.dT;                         // [hour]   Size of the time step which updates everything else
    dR                  = other.dR;                         // Depth of the infected layer (for spawning)
    nSamp               = other.nSamp;                      // Number of samples to save
    dGrid               = other.dGrid;                      // Spacing of grid (in units of critical radius R)
    dShell              = other.dShell;                     // [µm]     Thickness of agent layer (for biomass)

    n_0                 = other.n_0;                        // [1/µm^3] Initial concentration of nutrient
    g_max               = other.g_max;                      // [1/hour] Maximal growth rate for the cells
    R                   = other.R;                          // [µm]     The length scale for division (Typical volume 1.3 µm^3)
    k                   = other.k;                          // [N*m]    Parameter for repulsive potential
    K                   = other.K;                          //          Michaels-Menten FACTOR for Monod growth
    gamma               = other.gamma;                      //          Probability to infect cell
    alpha               = other.alpha;                      //          Probability for phage to go lysogenic
    beta                = other.beta;                       //          Multiplication factor phage
    delta               = other.delta;                      // [1/hour] Rate of phage decay
    epsilon             = other.epsilon;                    //          Probability for offspring to turn resistant
    r                   = other.r;                          //          Constant used in the time-delay mechanism
    T_i                 = other.T_i;                        // [hours]  Time when the phage infections begins (less than 0 disables phage infection)
    L_box               = other.L_box;                      // [µm]     Length of boundary condition box

    numB                = other.numB;                       // Current tally of suceptible cells
    numL                = other.numL;                       // Current tally of lysogenic cells
    numI                = other.numI;                       // Current tally of infected cells
    numD                = other.numD;                       // Current tally of dead cells
    numR                = other.numR;                       // Current tally of resistant cells

    eta                 = other.eta;                        // Amount of division noise (width of gaussian)
    nu                  = other.nu;                         // Amount of displacement noise during division. x = x0 + rand(nu) etc.

    D_B                 = other.D_B;                        // [µm^2/hour] Diffusion constant for the cells
    D_P                 = other.D_P;                        // [µm^2/hour] Diffusion constant for the phage
    D_n                 = other.D_n;                        // [µm^2/hour] Diffusion constant for the nutrient

    h_agar              = other.h_agar;                     // [µm]    Total height of the soft agar layers
    h_cell              = other.h_cell;                     // [µm]    Distance of the cell colony to the hard agar

    p                   = other.p;                          // Helper variable to control phage spawning

    debug               = other.debug;                      // The amount of information to print to terminal
    debugBool           = other.debugBool;
    exit                = other.exit;                       // Boolean to control early exit
    Time                = other.Time;                       // Counter for how many timesteps have passed
    maxStep             = other.maxStep;                    // Maximum time-step used in adaptive algorithm
    sigma               = other.sigma;                      // Maximum number of standard deviations to move in diffusions term
    posibleStep         = other.posibleStep;                // Largest possible timestep at current time

    nutrientField       = other.nutrientField;              // Boolean to toggle between nutrient field and just nutrients

    exportAny           = other.exportAny;                  //
    exportCellData      = other.exportCellData;             //
    exportColonySize    = other.exportColonySize;           // Booleans to control the export output
    exportPhageData     = other.exportPhageData;            //
    exportNutrient      = other.exportNutrient;             //

    singleInfectedCell  = other.singleInfectedCell;         // Boolean to enable phages via a single infected cell
    planarPhageInvasion = other.planarPhageInvasion;        // Boolean to enable phages invading in a plane
    uniformPhageInvasion= other.uniformPhageInvasion;       // Boolean to enable phages invading uniformly from entire space
    manyInfectedCells   = other.manyInfectedCells;          // Boolean to enable surface infection and uniformly from entire spacefalse
    forcedThickness     = other.forcedThickness;            // "Boolean" to force an infected surface of a given number

    invasionType        = other.invasionType;               // Number to remember the invasion type

    reflectingBoundary  = other.reflectingBoundary;         // Boolean to enable the simple box boundary condition with reflective sides
    absorbingBoundary   = other.absorbingBoundary;          // Boolean to enable the simple box boundary condition with absorbing sides
    experimentBoundary  = other.experimentBoundary;         // Boolean to enable the advanced petri dish boundary conditions

    phageType           = other.phageType;                  // Contains the type of phage chosen (sets parameters according to De Paepe)
    rngSeed             = other.rngSeed;                    // The seed for the random number generator

    // Initilize the benchmarking vector
    for (int j = 0; j < 20; j ++) {
        benchMark[j] = other.benchMark[j];
    }

    res = other.res;

    // Copy random number generator
    rng = other.rng;

    // Copy cells, phages, nutrients and nutrient grid.
    cells               = other.cells;
    biomass             = other.biomass;
    phages              = other.phages;
    center              = other.center;
    nutrient            = other.nutrient;
    nutrient_grid       = other.nutrient_grid;
    CN                  = other.CN;
    B                   = other.B;
}


// Controls the evaluation of the simulation
int Simulation::Run(double T) {

    // Check for invalid runs
    if (exit) { return 1; }
    if (T < dT) { return 1; }

    // Get start time
    time_t  tic;
    time(&tic);

    // Things to run only when simulation is initialized
    if (Time == 0) {

        // Initilize the simulation matrices
        Initialize();

    }

    // Check if export has been enabled, and if so, generate a path
    if (exportAny) {
        path = GeneratePath();

        // Write the reproducable command to log.txt
        WriteLog(T);
    }


    // Determine the number of timesteps between samplings
    int nStepsPerSample = (int)round(1/(nSamp*dt));

    // Determine the number of timesteps between cell updates
    int nStepsPerCellUpdate = (int)round(dT_c/dt);

    // Determine the number of samples to take during this run
    int nSamplings = nSamp*T;

    // Run steps which should not be sampled
    if ((Time < (int)round(T_i/dt)) or (nSamplings == 0)) {

        // If no samplings are to be made. Compute the number of time-steps to take:
        if (nSamplings == 0) {

            nStepsPerSample = (int)round(T/dt);

        } else {    // If simulation should run up until T_i, and then start sampling:

            // Compute number of steps between now and T_i
            nStepsPerSample = (int)round(T_i/dt) - Time - 1;
        }

        int t = 0;
        while (t < nStepsPerSample) {

            // Check for exit flag
            if (exit) { break; }

            // Update cells
            CellUpdate();

            // Use adaptive timesteps to update phages
            int tp = 0;
            while (tp < nStepsPerCellUpdate) {

                 // Adaptive timesteps
                int step = posibleStep;

                // Ensure the step is less than maxStep
                if (step > maxStep) {
                    step = maxStep;
                }

                // Ensure the step fits with next saving point
                if (step > nStepsPerCellUpdate-tp) {
                    step = nStepsPerCellUpdate-tp;
                }

                // Ensure the step is positive
                if (step < 1) {
                    step = 1;
                }

                dT = step*dt;
                PhageUpdate();

                // Increase Time counter and t if steps were stepped (t and Time are measured in units of dtF)
                Time   += step;
                tp     += step;
                t      += step;
                nSkips += step-1;
            }

            // if (numB == 0) {
            //     cerr << "\t>>No more uninfected cells! Exiting..<<" << endl;
            //     f_log << ">>No more uninfected cells! Exiting..<<" << endl;
            //     exit = true;
            // }

            // if (phages.size() == 0) {
            //     cerr << "\t>>No more phages! Exiting..<<" << endl;
            //     f_log << ">>No more phages! Exiting..<<" << endl;
            //     exit = true;
            // }

            if ((numB == 0) and (numI == 0) and (numL == 0)) {
                cerr << "\t>>All cells are dead! Exiting..<<" << endl;
                f_log << ">>All cells are dead! Exiting..<<" << endl;
                exit = true;
            }
        }

        // Reset the sampling frequency
        nStepsPerSample = (int)round(1/(nSamp*dt));
    }

    // Export the start configuration
    if (firstRun) {

        if ( exportCellData          and (debug > 0)) {cout << "\tExporting Cell Position Data" << endl;}
        if ( exportColonySize        and (debug > 0)) {cout << "\tExporting Colony Size" << endl;}
        if ( exportPhageData         and (debug > 0)) {cout << "\tExporting Phage Position Data" << endl;}
        if ( exportNutrient          and (debug > 0)) {cout << "\tExporting Nutrient" << endl;}

        // Export data
        if (exportAny) ExportData(Time*dt);

    }

    // Run the time evolution
    time_t timer;
    time(&timer);

    // Loop over samplings
    for (int n = 0; n < nSamplings; n++) {

        // Determine the number of samplings during run is compatible with timestep
        // cout << nSamplings << endl;
        // cout << nStepsPerSample << endl;
        // cout << dt << endl;
        // cout << T << endl;
        // cout << fabs(nSamplings*nStepsPerSample*dt-T) << endl;
        if (fabs(nSamplings*nStepsPerSample*dt-T) > dt/2) {
            cerr << "\t>>Time-step too large for sampling frequency! Exiting..<<" << endl;
            f_log << ">>Time-step too large for sampling frequency! Exiting..<<" << endl;
            exit = true;
        }

        // Check for exit flag
        if (exit) { break; }

        int t = 0;
        while (t < nStepsPerSample) {

            // Check for exit flag
            if (exit) { break; }

            // Update cells
            CellUpdate();

            // Use adaptive timesteps to update phages
            int tp = 0;
            while (tp < nStepsPerCellUpdate) {

                 // Adaptive timesteps
                int step = posibleStep;

                // Ensure the step is less than maxStep
                if (step > maxStep) {
                    step = maxStep;
                }

                // Ensure the step fits with next saving point
                if (step > nStepsPerCellUpdate-tp) {
                    step = nStepsPerCellUpdate-tp;
                }

                // Ensure the step is positive
                if (step < 1) {
                    step = 1;
                }

                dT = step*dt;
                PhageUpdate();

                // Increase Time counter and t if steps were stepped (t and Time are measured in units of dtF)
                Time   += step;
                tp     += step;
                t      += step;
                nSkips += step-1;
            }

            // if (numB == 0) {
            //     cerr << "\t>>No more uninfected cells! Exiting..<<" << endl;
            //     f_log << ">>No more uninfected cells! Exiting..<<" << endl;
            //     exit = true;
            // }

            // if (phages.size() == 0) {
            //     cerr << "\t>>No more phages! Exiting..<<" << endl;
            //     f_log << ">>No more phages! Exiting..<<" << endl;
            //     exit = true;
            // }

            if ((numB == 0) and (numI == 0) and (numL == 0)) {
                cerr << "\t>>All cells are dead! Exiting..<<" << endl;
                f_log << ">>All cells are dead! Exiting..<<" << endl;
                exit = true;
            }
        }


        // Export the data
        time(&timer);
        if (exportAny) ExportData(Time*dt);
        benchMark[6] += difftime(time(NULL),timer);

        // Show progress bar
        if ((n > 0) and (debug > 0)) {
            cout << "\t[";
            int pos = 60 * static_cast<float>(n)/static_cast<float>(nSamplings);;
            for (int i = 0; i < 60; ++i) {
                if (i <= pos) cout << ".";
                else cout << " ";
            }
            cout << "] " << "\r";
            cout.flush();
        }

        // Store the state to file
        time(&timer);
        if (exportAny) {
            SaveState();
        }
        benchMark[6] += difftime(time(NULL),timer);
    }

    // Get stop time
    time_t  toc;
    time(&toc);

    // Calculate time difference
    float seconds = difftime(toc, tic);
    float hours   = floor(seconds/3600);
    float minutes = floor(seconds/60);
    minutes -= hours*60;
    seconds -= minutes*60 + hours*3600;

    if (debug > 0) {
        cout << endl;
        cout << "\tSimulation complete after ";
        if (hours > 0.0)   cout << hours   << " hours and ";
        if (minutes > 0.0) cout << minutes << " minutes and ";
        cout  << seconds << " seconds." << endl;
    }

    // Report benchmarking results
    if (debug > 0) {
        cout << endl;
        cout << "\tBenchmarking results:" << endl;
        cout << "\t"       << setw(3) << difftime(toc, tic) << " s of total time" << endl;
        cout << "\t"       << setw(3) << benchMark[0]  << " s spent on computing movement" << endl;
        cout << "\t\t"     << setw(3) << benchMark[7]  << " s spent on computing cell bursting" << endl;
        cout << "\t\t"     << setw(3) << benchMark[8]  << " s spent on computing cell movement" << endl;
        cout << "\t\t"     << setw(3) << benchMark[9]  << " s spent on computing phage movement" << endl;
        cout << "\t\t\t"   << setw(3) << benchMark[12] << " s spent on phage-cell overlap and diffusion" << endl;
        cout << "\t\t\t\t" << setw(3) << benchMark[18] << " s spent on phage-cell overlap" << endl;
        cout << "\t\t\t\t" << setw(3) << benchMark[19] << " s spent on phage diffusion" << endl;
        cout << "\t\t\t"   << setw(3) << benchMark[13] << " s spent on phage-boundary reflections" << endl;
        cout << "\t"       << setw(3) << benchMark[1]  << " s spent on executing movement" << endl;
        cout << "\t"       << setw(3) << benchMark[2]  << " s spent on cell growth" << endl;
        cout << "\t"       << setw(3) << benchMark[3]  << " s spent on nutrient" << endl;
        cout << "\t"       << setw(3) << benchMark[4]  << " s spent on infections" << endl;
        cout << "\t"       << setw(3) << benchMark[5]  << " s spent on CN, B and P grids" << endl;
        cout << "\t\t"     << setw(3) << benchMark[10] << " s spent on armadillo data structures" << endl;
        cout << "\t\t\t"   << setw(3) << benchMark[14] << " s spent on clearing arrays" << endl;
        cout << "\t\t\t"   << setw(3) << benchMark[15] << " s spent on updating arrays" << endl;
        cout << "\t\t"     << setw(3) << benchMark[11] << " s spent on std vector data structures" << endl;
        cout << "\t\t\t"   << setw(3) << benchMark[16] << " s spent on clearing arrays" << endl;
        cout << "\t\t\t"   << setw(3) << benchMark[17] << " s spent on updating arrays" << endl;
        cout << "\t"       << setw(3) << benchMark[6]  << " s spent on filesaving" << endl;

        cout << "\t----------------------------------------------------"<< endl << endl << endl;
    }

    // Write sucess to log
    if (exit) {
        f_log << "Adaptive algorithm reduced time-steps by: " << nSkips << endl;
        f_log << ">>Simulation completed with exit flag<<" << endl;
    } else {
        f_log << "Adaptive algorithm reduced time-steps by: " << nSkips << endl;
        f_log << ">>Simulation completed without exit flag<<" << endl;
    }

    if (exit) {
        return 1;
    } else {
        return 0;
    }
}


// Initialize the simulation
void Simulation::Initialize() {

    // Set the random number generator seed
    if (rngSeed >= 0.0) {
        rng.seed( rngSeed );
    } else {
        static std::random_device rd;
        rng.seed(rd());
    }

    // Initilize unit nutrient throughout the simulation
    if (debug > 1) {
        cout << endl;
        cout << "\tInitlizing the nutrient grid, cell density grid, phage density grid, clever neighbour grid, and laplace operator" << endl;
    }

    // Resize armadillo data types
    if (nutrientField) {
        nutrient_grid.set_size(2*res+1,2*res+1,2*res+1);
    }

      B.set_size(2*res+1,2*res+1,2*res+1);
    lap.set_size(2*res+1,2*res+1);

    // Fill nutrient and cell density
    if (nutrientField) nutrient_grid.fill(n_0*pow(dGrid*R,3));  // Nutrient per unit cell
    else nutrient = n_0*pow(dGrid*R,3);                         // Average nutrient per unit cell
    B.zeros();

    // Compute the proper value for the Michaels-Menten
    K *= n_0*pow(dGrid*R,3);

    // Reinterpret dGrid
    dGrid = dGrid*R;

    // Adjust the maxStep
    // maxStep = min(maxStep,(int)floor(pow(dGrid*R,2)/(pow(1,2)*6*D_P*dT)));

    // Boundary conditions for laplace operator
    lap(0,0)            = -1;
    lap(0,1)            =  1;
    lap(2*res,2*res)    = -1;
    lap(2*res,2*res-1)  =  1;


    // Fill the laplace operator, and resize the clever neighbour grid
    CN.resize(2*res+1);
    for (int x = 0; x < 2*res+1; x++ ) {

        if ((x > 0) and (x < 2*res)) {
            lap(x,x)    = -2;
            lap(x,x+1)  =  1;
            lap(x,x-1)  =  1;
        }

        CN[x].resize(2*res+1);
        for (int y = 0; y < 2*res+1; y++ ) {

            CN[x][y].resize(2*res+1);
        }
    }

    // Initilize the cells at (0,0,0) + random offset, with radius r = 0.677+random offset, and in the 0 state
    if (debug > 1) {cout << "\tInitlizing the cells" << endl;}

    cells.reserve(N_max);
    cells.resize(N_0);
    for (int n = 0; n < N_0; n++) {

        // Generate location
        double x = (2*rand(rng)-1)*N_0*R;
        double y = (2*rand(rng)-1)*N_0*R;
        double z = (2*rand(rng)-1)*N_0*R;

        // Generate radius
        double rad = 0.0;
        if (n == 0) {
            rad = 0.677;
        } else {
            rad = 0.677*(1 + (Rand(0.2)-0.1));
        }


        // Detect which grid point cell belongs to
        int i = round(x / (dGrid)) + res;
        int j = round(y / (dGrid)) + res;
        int k = round(z / (dGrid)) + res;

        // Add cell to system
        cells[n] = vector<double> { x, y, z, rad, 0, (double)i, (double)j, (double)k, 0};

        // Add information to clever neighbour grid
        CN[i][j][k].push_back(n);

        // Update tally
        numB++;

    }

    // Add center of colony
    center = vector<double> { 0, 0, 0, -1 };
    r_max = 0;
    I_max = -1;

    // Do some checks for input parameters
    if ( (not exit) and (nutrientField) and (2*D_n *dT_c/pow(dGrid,2) > 1) ) {
        cerr << "\t>>Unstable resolution sizes [ 2*D_n dT/(dGrid)^2 = " << 2*D_n*dT_c/pow(dGrid,2) << " > 1 ]! Exiting..<<" << endl;
        f_log << ">>Unstable resolution sizes [ 2*D_n dT/(dGrid)^2 = " << 2*D_n*dT_c/pow(dGrid,2) << " > 1 ]! Exiting..<<" << endl;
        exit = true;
    }
}


// Updates the cells
void Simulation::CellUpdate() {

    // Check for exit flag
    if (exit) { return; }

    // Time objects for benchmarking
    time_t timer;

    // Reset tally of uninfected
    numB = 0;

    if (debugBool) {deb(11);}
    time(&timer);
    // Increase the state of infections, and spawn new phages.
    bool bursting = false;
    int N = cells.size();
    for (int n = N-1; n >= 0 ; n--) {

        if (GrowInfection(n)) {

            // Remove bursted cell
            cells.erase(cells.begin() + n);

            // Set bursting variable to true
            bursting = true;

            // Update tally
            numI--;
            numD++;
        }
    }

    benchMark[7] += difftime(time(NULL),timer);
    benchMark[0] += difftime(time(NULL),timer);

    // Check if all cells are bursted
    if (cells.size() == 0) {
        cerr << "\t>>All cells have bursted! Exiting..<<" << endl;
        f_log << ">>All cells have bursted! Exiting..<<" << endl;
        exit = true;
        return;
    }

    // Check for exit flag
    if (exit) { return; }


    if (debugBool) {deb(12);}
    time(&timer);
    // Grow the cells
    N = cells.size();
    for (int n = 0; n < N; n++) {

        if (cells[n][4] > 0) {              // Lytic cells
            continue;
        } else if (cells[n][4] == 0) {      // Uninfected cells
            numB++;
            GrowCell(n);
        } else if (cells[n][4] == -1) {     // Lysogenic cells
            GrowCell(n);
        } else if (cells[n][4] == -2) {     // Dead cells
            continue;
        } else if (cells[n][4] == -3) {     // Resistant cells
            GrowCell(n);
        }

    }

    if (debugBool) {deb(13);}
    // Grow the biomass
    if (not biomass.empty()) {
        if (GrowBiomass()) {
            // If biomass absorbed a cell, update CN grid
            bursting = true;
        }

        // Update the number of uninfected cells
        numB += (int)floor(pow(biomass[3],3)/pow(0.677,3));
    }

    benchMark[2] += difftime(time(NULL),timer);

    if (debugBool) {deb(14);}
    // Update the CN grid if cells have been removed
    if (bursting) {
        UpdateNearestNeighbourGrid();
    }

    if (debugBool) {deb(15);}
    if (cells.size() >= 100) {
        if (omp_get_num_threads() == 1) {
            omp_set_num_threads(omp_get_max_threads());
        }
    }

    if (debugBool) {deb(16);}
    time(&timer);
    // Run the over cells to get movement for the cells
    N = cells.size();
    double** C_movement = new double*[N];
    #pragma omp parallel for
    for (int n = 0; n < N; n++) {
        C_movement[n] = new double[3];
        CellMovement(n,C_movement[n]);
    }
    benchMark[8] += difftime(time(NULL),timer);

    if (debugBool) {deb(17);}
    // Center the biomass
    if (not biomass.empty()) {
        int counter = 0;
        double ucenter[3] = {0.0, 0.0, 0.0}; // Center of uninfected cells
        for (int n = 0; n < N; n++) {
            if (cells[n][4] == 0) {
                for (int j = 0; j < 3; j ++) {
                    ucenter[j] += cells[n][j];
                    counter++;
                }
            }
        }

        // Move the biomass
        for (int j = 0; j < 3; j ++) {
            biomass[j] = 0.25*ucenter[j]/counter + (1-0.25)*biomass[j];
        }

        // Detect which grid point center belongs to
        biomass[5] = round( biomass[0] / (dGrid) ) + res;
        biomass[6] = round( biomass[1] / (dGrid) ) + res;
        biomass[7] = round( biomass[2] / (dGrid) ) + res;
    }


    if (debugBool) {deb(18);}
    time(&timer);
    // Update the nutrient
    if (K > 0.0) {
        if (nutrientField) {
            NutrientGridUpdate();
            if (accu(nutrient_grid)/(pow(2*res+1,3)*pow(dGrid,3)) < 1e-9) {
                cerr << "\t>>Nutrients are depleted! (Time = " << Time*dt << "); Exiting..<<" << endl;
                f_log << ">>Nutrients are depleted! (Time = " << Time*dt << "); Exiting..<<" << endl;
                exit = true;
            }
        } else {
            nutrient -= g_max*nutrient/(nutrient+K)*(numB+numL)/pow(2*res+1,3)*dT;
            nutrient = fmax(0.0,nutrient); // Ensure that negative numbers does not occur.
            if (nutrient/pow(dGrid,3) < 1e-9) {
                cerr << "\t>>Nutrients are depleted! (Time = " << Time*dt << "); Exiting..<<" << endl;
                f_log << ">>Nutrients are depleted! (Time = " << Time*dt << "); Exiting..<<" << endl;
                exit = true;
            }
        }
    }
    benchMark[3] += difftime(time(NULL),timer);


    if (debugBool) {deb(19);}
    time(&timer);
    // Reset the cell density grid and phage grid
    if (nutrientField) {
        B.fill(0);
    }
    benchMark[14] += difftime(time(NULL),timer);
    benchMark[10] += difftime(time(NULL),timer);
    benchMark[5]  += difftime(time(NULL),timer);

    if (debugBool) {deb(20);}
    // Run the over cells to apply the movement step
    for (int n = 0; n < N; n++) {

        time(&timer);
        for (int j = 0; j < 3; j++) {
            cells[n][j] += C_movement[n][j];
        }
        benchMark[1] += difftime(time(NULL),timer);

        // Detect which grid point cell belongs to now
        int i = (int)round(cells[n][0] / (dGrid)) + res;
        int j = (int)round(cells[n][1] / (dGrid)) + res;
        int k = (int)round(cells[n][2] / (dGrid)) + res;

        // Check if outer bondary is reached
        if ( (i < 0) or (i >= 2*res+1) or (j < 0) or (j >= 2*res+1) or (k < 0) or (k >= 2*res+1) ) {
            cerr << "\t>>Colony extends simulation grid size! Exiting...<<" << endl;
            f_log << ">>Colony extends simulation grid size! Exiting...<<" << endl;
            exit = true;
            return;
        }

        // Get the old gridpoint
        int io = (int)cells[n][5];
        int jo = (int)cells[n][6];
        int ko = (int)cells[n][7];

        // Check if cell moved to new point
        if ( (not ( io == i)) or (not ( jo == j)) or (not ( ko == k))) {

            // Erase from old grid point
            CN[io][jo][ko].erase(std::remove(CN[io][jo][ko].begin(), CN[io][jo][ko].end(), n), CN[io][jo][ko].end());

            // Update location information
            cells[n][5] = (double)i;
            cells[n][6] = (double)j;
            cells[n][7] = (double)k;

            // Add cell to new CN grid point
            CN[i][j][k].push_back(n);

        }

        // Add cell to cell density grid
        if (nutrientField) {
            if ((cells[n][4] == 0) or (cells[n][4] == -1)) {
                B(i,j,k)++;
            }
        }
        benchMark[15] += difftime(time(NULL),timer);
        benchMark[10] += difftime(time(NULL),timer);
        benchMark[17] += difftime(time(NULL),timer);
        benchMark[11] += difftime(time(NULL),timer);
    }


    if (debugBool) {deb(21);}
    // Replaces the most central bacteria with a central biomass if the colony is too large
    if (not exit) {
        if (biomass.empty()) {
            CreateCentralBiomass();
        }
    }

    if (debugBool) {deb(22);}
    // Update the colony extent
    ComputeColonyExtent();

    // Clean up
    if (debugBool) {deb(23);}
    time(&timer);
    for (int n = 0; n < N; n++) {
        delete[] C_movement[n];
    }
    delete[] C_movement;

    if (debugBool) {deb(24);}
}





// Updates the phages
void Simulation::PhageUpdate() {

    // Check for exit flag
    if (exit) { return; }

    // Time objects for benchmarking
    time_t timer;

    if (debugBool) {deb(1);}
    time(&timer);
    // Run the over phages to infect the cells
    int M = phages.size();
    for (int m = M-1; m >= 0; m--) {

        // Remove the phage if it sucessfully infects the cell
        if (PhageInfection(m) == 1) {

            phages.erase(phages.begin() + m);

        }

    }
    benchMark[4] += difftime(time(NULL),timer);

    if (debugBool) {deb(2);}
    time(&timer);
    // Run the over phages to get movement for the phage
    M = phages.size();
    double** P_movement = new double*[M];
    #pragma omp parallel for
    for (int m = 0; m < M; m++) {

        P_movement[m] = new double[3];
        PhageMovement(m,P_movement[m]);

    }
    benchMark[9]  += difftime(time(NULL),timer);
    benchMark[0]  += difftime(time(NULL),timer);


    if (debugBool) {deb(3);}
    time(&timer);
    // Run the over phages to apply the movement step
    for (int m = M-1; m >= 0; m--) {

        time(&timer);
        // Kill phages according to the decay rate
        if (rand(rng) < delta*dT) {

            phages.erase(phages.begin() + m);
            continue;

        }

        // Apply the movement step
        for (int j = 0; j < 3; j++) {
            phages[m][j] += P_movement[m][j];
        }

        benchMark[1] += difftime(time(NULL),timer);

        time(&timer);
        // Apply boundary conditions to the phage
        ApplyBoundaryConditions(m);

        benchMark[13] += difftime(time(NULL),timer);
        benchMark[9]  += difftime(time(NULL),timer);
        benchMark[0]  += difftime(time(NULL),timer);
    }

    if (debugBool) {deb(4);}

    // Spawns phages, if the infection time has been passed
    if (not exit) {
        SpawnPhages();
    }

    // Clean up
    for (int m = 0; m < M; m++) {
        delete[] P_movement[m];
    }
    delete[] P_movement;

    if (debugBool) {deb(5);}
}


// Simulation functions /////////////////////////////////////////////////////////////////
// Get movement vector for cell I
void Simulation::CellMovement(int I, double *B) {

    // Prepare the output vector
    std::fill(B,B+3,0);

    // Get the cell gridpoint
    int i = cells[I][5];
    int j = cells[I][6];
    int k = cells[I][7];

    // Find nearest neighbours
    vector<int> NN = NearestNeighbours(i,j,k);
    for (int i = 0; i < NN.size(); i++) {

        // Neighbour to test overlap with
        int n = NN[i];

        // Skip the current cell ( cell I )
        if ( n == I ) continue;

        // Calculate distance between cells
        double d = 0.0;
        for (int j = 0; j < 3; j++) {
            d += pow(cells[n][j]-cells[I][j],2);
        }
        d = sqrt(d);                        // Distance between center of mass

        // Check for overlap
        if ( d < cells[n][3]+cells[I][3] ) {

            // Calculate the displacement of the cell
            // See notes p. 12
            double F = - PotentialGradient(d,cells[n][3],cells[I][3]) * dT_c;

            // Use the difference in position vectors to determine the direction of the force
            for (int j = 0; j < 3; j++) {
                B[j] += F * ( cells[n][j] - cells[I][j] ) / d;
            }
        }
    }

    // Check for overlap with central biomass
    if (not biomass.empty()) {

        // Calculate distance between cell and biomass
        double d = 0.0;
        for (int j = 0; j < 3; j++) {
            d += pow(cells[I][j]-biomass[j],2);
        }
        d = sqrt(d);                        // Distance between center of masses


        // Check for overlap
        if ( d < cells[I][3]+biomass[3]/pow(0.64,1/3) ) {

            // The overlap potential assumes overlap between cells of similar sizes, and
            // since the biomass is potentially huge compared to the cells, some modifications
            // must be made. I wish to model the dynamics of a single cell of radius R pushing the
            // bacteria back, rather than the biomass as one large cell.

            // Calculate the displacement of the cell
            // assert(d-biomass[3]/pow(0.64,1/3)+R > 0);
            // double F = - 100*PotentialGradient(d-biomass[3]/pow(0.64,1/3)+R,R,cells[I][3]) * dT;
            // double F = - floor(pow(biomass[3],3)/pow(0.677,3))*PotentialGradient(d,biomass[3]/pow(0.64,1/3),cells[I][3]) * dT;
            double F = - 10*PotentialGradient(d,biomass[3]/pow(0.64,1/3),cells[I][3]) * dT_c;

            // Use the difference in position vectors to determine the direction of the force
            for (int j = 0; j < 3; j++) {
                B[j] += F * ( cells[I][j]-biomass[j] ) / d;
            }

        }

    }

    // Check if diffusion is enabled
    if (D_B > 0) {

        // Calculate the change due to the diffusion
        for (int j = 0; j < 3; j++) {
            B[j] += sqrt(2*D_B*dT_c)*randn(rng);
        }
    }
}


// Get movement vector for phage I
void Simulation::PhageMovement(int I, double *P) {

    // Prepare the output vector
    std::fill(P,P+3,0);

    time_t timer;
    time(&timer);

    // Detect which grid point cell belongs to
    int i = round(phages[I][0] / (dGrid)) + res;
    int j = round(phages[I][1] / (dGrid)) + res;
    int k = round(phages[I][2] / (dGrid)) + res;

    // Find nearest neighbours
    vector<int> NN = NearestNeighbours(i,j,k);

    for (int i = 0; i < NN.size(); i++) {

        // Neighbour to test overlap with
        int n = NN[i];

        // Calculate distance between phage and cell
        double d = 0.0;
        for (int j = 0; j < 3; j++) {
            d += pow(cells[n][j]-phages[I][j],2);
        }
        d = sqrt(d);                        // Distance between center of mass


        // Check for overlap
        if ( d < cells[n][3] ) {

            // Calculate the displacement of the phage (using same potential as the cells)
            // See notes p. 12
            double F = - PotentialGradient(d,cells[n][3],0) * dT;

            // Use the difference in position vectors to determine the direction of the force
            for (int j = 0; j < 3; j++) {
                P[j] += F * ( phages[I][j]-cells[n][j] ) / d;
            }
        }

        // Set the adaptive timescale
        // Phage is "close" to colony (ie. there are cells within NN)

        // Compute the number of time-steps it takes to move the distance from phage to current cell
        // int s = max(1,(int)ceil(1*pow(d-cells[n][3],2)/(4*M_PIl*D_P*dt))); // t ~ L^2/D
        int s = max(1,(int)floor(pow((d-cells[n][3])/sigma,2)/(2*D_P*dt)));    // x ~ sqrt(2*D*dT)

        if (s < posibleStep) {posibleStep = s;}
    }


    // Set the adaptive timescale
    // Phage is "far" from colony (ie. there are NO cells within NN)

    // Compute distance from phage to center of colony
    double d = sqrt(pow(phages[I][0]-center[0],2)+pow(phages[I][1]-center[1],2)+pow(phages[I][2]-center[2],2));

    // Compute distance from phage to edge of colony
    d = d - r_max;

    // If phage is outside the colony extent
    if ( d > 0 ) {
        // Compute the number of time-steps it takes to move the distance from phage to edge of colony
        // int s = max(1,(int)ceil(1*pow(d,2)/(4*M_PIl*D_P*dt))); // t ~ L^2/D
        int s = max(1,(int)floor(pow(d/sigma,2)/(2*D_P*dt)));    // x ~ sqrt(2*D*dT)
        if (s < posibleStep) {posibleStep = s;}

    } else {    // Step based on NN grid size
        // Compute the number of time-steps it takes to move to next grid-point
        // int s = max(1,(int)ceil(1*pow(dGrid,2)/(4*M_PIl*D_P*dt))); // t ~ L^2/D
        int s = max(1,(int)floor(pow(dGrid/sigma,2)/(2*D_P*dt)));    // x ~ sqrt(2*D*dT)
        if (s < posibleStep) {posibleStep = s;}
    }


    // Check if diffusion is enabled
    if (D_P > 0) {

        // Calculate the change due to the diffusion
        for (int j = 0; j < 3; j++) {
            double dx = randn(rng);
            if (dx > 0) {
                dx = min(sigma,dx);
            } else {
                dx = max(-sigma,dx);
            }
            P[j] += sqrt(2*D_P*dT)*dx;
        }
    }

    benchMark[18] += difftime(time(NULL),timer);
}


// Check for infection oppotunity for phage I
int Simulation::PhageInfection(int I) {

    // Get the grid coordinate for the phage
    int i = round(phages[I][0] / (dGrid)) + res;
    int j = round(phages[I][1] / (dGrid)) + res;
    int k = round(phages[I][2] / (dGrid)) + res;

    // Ensure coordinate is within CN area
    if ( (i < 0) or (i > 2*res+1) or (j < 0) or (j > 2*res+1) or (k < 0) or (k > 2*res+1) ) {
        return 0;
    }

    // Find nearest neighbours
    vector<int> NN = NearestNeighbours(i,j,k);

    for (int i = 0; i < NN.size(); i++) {

        // Neighbour to test overlap with
        int n = NN[i];

        // Calculate distance between phage and cell
        double d = 0.0;
        for (int l = 0; l < 3; l++) {
            d += pow(cells[n][l]-phages[I][l],2);
        }
        d = sqrt(d);

        // Check if phage is inside cell
        if (d < cells[n][3]) {


            // Pass a probability check
            if (rand(rng) <= gamma*dT) {

                // Count the infection
                cells[n][8]++;

                // Increment the phage state if cell is uninfected
                if (cells[n][4] == 0) {

                    // Choose infection pathway
                    if (rand(rng) <= alpha) {  // Lysogenic
                        cells[n][4]--;
                        numL++;
                        // numB--;
                    } else {                // Lytic
                        cells[n][4]++;
                        numI++;
                        // numB--;
                    }

                    // return "success" flag
                    return 1;
                }


                // Absorb the phage if cell is lytic or lysogenic
                if ((cells[n][4] == -1) or (cells[n][4] >= 1)) {

                    // return "success" flag
                    return 1;
                }

                // Ensure phage is not absorbed by dead cells
                if (cells[n][4] == -2) {

                    // return "failure" flag
                    return 0;
                }

                // Ensure phage is not absorbed by resistant cells
                if (cells[n][4] == -3) {

                    // return "failure" flag
                    return 0;
                }


            } else {

                // return "failure" flag
                return 0;
            }
        }
    }

    // If no succesfull infection have occured, return 0.
    return 0;
}


// Update radius for a single cell (and grow new cells)
void Simulation::GrowCell(int I) {

    // Extract the growthrate for the location
    double g = 0.0;
    if (nutrientField) g = GrowthRate( cells[I][0], cells[I][1], cells[I][2] );     // Get available nutrient in unit cell
    else if (K > 0.0) g = g_max * nutrient/(nutrient+K);
    else g = g_max;

    // If the cell is too far away from the edge, grow the biomass instead
    if (not biomass.empty()) {
        double r = sqrt(pow(cells[I][0]-center[0],2) + pow(cells[I][1]-center[1],2) + pow(cells[I][2]-center[2],2));
        if ((r_max - r) > dShell) {
            biomass[3] = pow(pow(biomass[3],3)+pow(dT_c*g/3.0*cells[I][3],3),1.0/3.0);
        }
        return;
    }

    // Assuming exponential growth
    cells[I][3] += dT_c * g/3.0 * ( cells[I][3] );


    // Divide cell
    if (cells[I][3] > R) {

        // Spawn new cell
        if (cells.size() < N_max) {

            // Division leads uneven distribution of volume
            //      V -> alpha*V + (1-alpha)*V
            // =>   alpha*(r)^3 -> (r')^3
            // =>   r1' = alpha*r0
            // and  r2' = (1-alpha)*r0

            // Draw alpha from normal distribution
            double alpha = -1.0;
            while ((alpha < 0) or (alpha > 1)) {
                alpha = RandN(0.5,eta);
            }

            // Make uneven distribution of volume
            double r1 = cells[I][3] * cbrt(alpha);
            double r2 = cells[I][3] * cbrt(1-alpha);

            // Choose random angle
            double theta =   M_PIl*rand(rng);
            double phi   = 2*M_PIl*rand(rng);

            // Generate endpoint cooordinates
            double x1 = cells[I][0] + r1*sin(theta)*cos(phi);
            double y1 = cells[I][1] + r1*sin(theta)*sin(phi);
            double z1 = cells[I][2] + r1*cos(theta);

            double x2 = cells[I][0] + r2*sin(M_PIl-theta)*cos(phi+M_PIl);
            double y2 = cells[I][1] + r2*sin(M_PIl-theta)*sin(phi+M_PIl);
            double z2 = cells[I][2] + r2*cos(M_PIl-theta);

            // Update position and radius of old cell
            cells[I][0] = x1;
            cells[I][1] = y1;
            cells[I][2] = z1;
            cells[I][3] = r1;

            // Get current gridpoint of old cell
            int io = (int)cells[I][5];
            int jo = (int)cells[I][6];
            int ko = (int)cells[I][7];

            // Compute new gridpoint of old cell
            int i = round(x1 / dGrid) + res;
            int j = round(y1 / dGrid) + res;
            int k = round(z1 / dGrid) + res;

            // Check if cell moved to new point
            if ( (not ( io == i)) or (not ( jo == j)) or (not ( ko == k))) {

                // Erase from old grid point
                CN[io][jo][ko].erase(std::remove(CN[io][jo][ko].begin(), CN[io][jo][ko].end(), I), CN[io][jo][ko].end());

                // Update location information
                cells[I][5] = (double)i;
                cells[I][6] = (double)j;
                cells[I][7] = (double)k;

                // Add cell to new CN grid point
                CN[i][j][k].push_back(I);

            }


            // Detect which grid point the new cell belongs to
            i = round(x2 / dGrid) + res;
            j = round(y2 / dGrid) + res;
            k = round(z2 / dGrid) + res;

            // Does the offspring mutate?
            double S = cells[I][4];
            if ((S == 0) and (rand(rng) < epsilon)) {
                S = -3;
            }

            // Add the new cell
            cells.push_back( vector<double> { x2, y2, z2, r2, S, (double)i, (double)j, (double)k, 0 } );

            // Add cell to CN grid
            CN[i][j][k].push_back(cells.size()-1);

            // Update tally
            if (S == -1) {
                numL++;
            } else if (S == -3) {
                numR++;
            } else {
                numB++;
            }

        } else if (not exit) {

            cerr << "\t>>Maximum number of cells reached! Exiting..<<" << endl;
            f_log << ">>Maximum number of cells reached! Exiting..<<" << endl;
            exit = true;
        }
    }
}


// Update radius for the biomass and spawn new cells
bool Simulation::GrowBiomass() {

    // Extract the growthrate for the location
    double g = 0.0;
    if (nutrientField) g = GrowthRate( biomass[0], biomass[1], biomass[2] );     // Get available nutrient in unit cell
    else if (K > 0.0) g = g_max * nutrient/(nutrient+K);
    else g = g_max;

    // Assuming exponential growth
    biomass[3] += dT_c * g/3.0 * ( biomass[3] );

    // Compute the maximal distance from central biomass to any infected bacteria
    double r_min = HUGE_VAL;
    for (int n = 0; n < cells.size(); n++) {

        double r = sqrt(pow(cells[n][0]-biomass[0],2) + pow(cells[n][1]-biomass[1],2) + pow(cells[n][2]-biomass[2],2));

        if (cells[n][4] > 0) {
            if (r < r_min) {
                r_min = r;
            }
        }
    }

    // Set r_min equal to colony extent if there are no infected bacteria
    if (numI == 0) {
        r_min = r_max;
    }



    // Check if this distance is too small and new cells should be created from biomass
    if (r_min - biomass[3]/pow(0.64,1/3) < dShell) {

        // Create new from biomass
        if (biomass[3] > R) {

            // Choose random angle and radius
            double theta =   M_PIl*rand(rng);
            double phi   = 2*M_PIl*rand(rng);
            double r     = biomass[3]/pow(0.64,1/3);

            double x = biomass[0] + r*sin(theta)*cos(phi);
            double y = biomass[1] + r*sin(theta)*sin(phi);
            double z = biomass[2] + r*cos(theta);

             // Detect which grid point cell belongs to
            int i = round(x / (dGrid)) + res;
            int j = round(y / (dGrid)) + res;
            int k = round(z / (dGrid)) + res;

            // Generate radius
            double rad = 0.677;

            // Add the new cell
            cells.push_back( vector<double> { x, y, z, rad, 0, (double)i, (double)j, (double)k } );

            // Add cell to CN grid
            CN[i][j][k].push_back(cells.size()-1);

            // Remove mass from the biomass
            biomass[3] = pow(pow(biomass[3],3)-pow(rad,3),1.0/3.0);

        // If there is too little biomass, re-interpret biomass a cell.
        } else {

            // Get location of biomass
            int i = round(biomass[0] / (dGrid)) + res;
            int j = round(biomass[1] / (dGrid)) + res;
            int k = round(biomass[2] / (dGrid)) + res;

            biomass[5] = (double)i;
            biomass[6] = (double)j;
            biomass[7] = (double)k;

            // Add biomass to cell list
            cells.push_back( biomass );

            // Update CN grid
            CN[i][j][k].push_back(cells.size()-1);

            // Delete biomass
            biomass.resize(0);

        }
    }

    return false;
}


// Update stage of infection, and create new phages.
bool Simulation::GrowInfection(int I) {

    // Extract the growthrate for the location
    double g = 0.0;
    if (nutrientField) g = GrowthRate( cells[I][0], cells[I][1], cells[I][2] );     // Get available nutrient in unit cell
    else if (K > 0.0) g = g_max * nutrient/(nutrient+K);
    else g = g_max;


    // Check if cell is infected and if infection stage should be increased
    if ((cells[I][4] > 0) and (rand(rng) <= r*g*dT_c)) {
        cells[I][4]++;
    }

    // Check if cell is ready to burst
    if (cells[I][4] == 11) {

        // Resize the phage array
        int M = phages.size();

        // Check if this would spawn too many phages
        if (M + beta > M_max) {

            cerr << "\t>>Maximum number of phages reached! Exiting..<<" << endl;
            f_log << ">>Maximum number of phages reached! Exiting..<<" << endl;
            exit = true;

        } else {

            phages.resize(M+beta);
            M_tot += beta;

            // Spawn new phages (Uniformly on the surface of the bursted cell)
            for (int b = 0; b < beta; b++) {

                // Choose random angle and radius
                double theta =       M_PIl*rand(rng);
                double phi   =     2*M_PIl*rand(rng);
                double r     = cells[I][3]*rand(rng)*0.90;

                double x = cells[I][0] + r*sin(theta)*cos(phi);
                double y = cells[I][1] + r*sin(theta)*sin(phi);
                double z = cells[I][2] + r*cos(theta);

                phages[M+b] = vector<double> { x, y, z } ;

            }

            // Set the the burst cell flag
            return true;
        }

    }

    // Cell does not burst
    return false;
}


// Returns the growthrate at the location
double Simulation::GrowthRate(double x, double y, double z) {

    // Detect which grid point cell belongs to
    int i = round( x / (dGrid) ) + res;
    int j = round( y / (dGrid) ) + res;
    int k = round( z / (dGrid) ) + res;

    // Check that cell is inside boundary
    if ( (i < 0) or (i > 2*res) or (j < 0) or (j > 2*res) or (k < 0) or (k > 2*res) ) {
        cout << "i = " << i << "; j = " << j << "; k = " << k << endl;
        cerr << "\t>>Cell has been pushed out of bounds! Crashing..<<" << endl;
        f_log << ">>Cell has been pushed out of bounds! Crashing..<<" << endl;
        exit = true;
    }

    // Assuming Monod growth law
    // g(n) = g_max * n / (n + K)

    // g_max is per minute
    return g_max*nutrient_grid(i,j,k)/(nutrient_grid(i,j,k) + K);
}


// Update the nutrient from diff. eqns.
void Simulation::NutrientGridUpdate() {

    // Create copy of nutrient to store the diffusion update
    cube nn = nutrient_grid;

    // Apply the growth rate step
    nutrient_grid -= g_max * (nutrient_grid/(nutrient_grid+K)) % B * dT_c;


    // Compute the X & Y diffusion
    for (int k = 0; k < 2*res+1; k++) {
        nutrient_grid.slice(k) += D_n*dT_c/pow(dGrid,2) * ( (lap*nn.slice(k).t()).t() + lap*nn.slice(k) );
    }

    // Compute the Z diffusion
    for (int i = 0; i < 2*res+1; i++) {

        mat Q = D_n*dT_c/pow(dGrid,2) * (lap*static_cast<mat>(nn.tube( span(i), span::all )).t()).t();


        for (int j = 0; j < Q.n_rows; j++) {

            for (int k = 0; k < Q.n_cols; k++) {
                nutrient_grid(i,j,k) += Q(j,k);
            }
        }
    }
}


// Applies the boundary conditions to the phages
void Simulation::ApplyBoundaryConditions(int m) {

    // Check X-direction
    if (phages[m][0] >  L_box/2) {
        if (absorbingBoundary) {phages.erase(phages.begin() + m); return;}
        phages[m][0] -= 2*(phages[m][0] - L_box/2);
    }
    if (phages[m][0] < -L_box/2) {
        if (absorbingBoundary) {phages.erase(phages.begin() + m); return;}
        phages[m][0] += 2*(-L_box/2 - phages[m][0]);
    }

    // Check Y-direction
    if (phages[m][1] >  L_box/2) {
        if (absorbingBoundary) {phages.erase(phages.begin() + m); return;}
        phages[m][1] -= 2*(phages[m][1] - L_box/2);
    }
    if (phages[m][1] < -L_box/2) {
        if (absorbingBoundary) {phages.erase(phages.begin() + m); return;}
        phages[m][1] += 2*(-L_box/2 - phages[m][1]);
    }

    // Check Z-direction
    if (reflectingBoundary or absorbingBoundary) {

        if (phages[m][2] >  L_box/2) {
            if (absorbingBoundary) {phages.erase(phages.begin() + m); return;}
            phages[m][2] -= 2*(phages[m][2] - L_box/2);
        }
        if (phages[m][2] < -L_box/2) {
            if (absorbingBoundary) {phages.erase(phages.begin() + m); return;}
            phages[m][2] += 2*(-L_box/2 - phages[m][2]);
        }

    } else if (experimentBoundary) {

        if (phages[m][2] > h_agar - h_cell) {   // Reflecting top boundary
            phages[m][2] -= 2*(phages[m][2] - h_agar + h_cell);
        } else if (phages[m][2] < -h_cell) { // Absorbing hard agar
            phages.erase(phages.begin() + m);
        }
    }
}


// Spanws M phages according to the spawning rules
void Simulation::SpawnPhages() {

    if (Time*dt <= T_i) {return;}

    // Check which spawning rule is used
    if (planarPhageInvasion) {

        // Compute number of phages to spawn
        p += FirstPassageTimeDistribution() * dT * (int)(P_0*pow(L_box,2));

        // Spawn phages
        while (p > 1) {

            // Random entry along XY plane (Uniformly)
            double x = (rand(rng)-0.5)*L_box;
            double y = (rand(rng)-0.5)*L_box;

            // Spawn the new phage
            phages.push_back( vector<double> {x, y, (res+0.5)*dGrid });
            M_tot++;

            // Reduce the number of phages to spawn.
            p -= 1.0;

        }

    }
    if ((singleInfectedCell) and (P_0 > 0)) {

        // Update the colony extent
        ComputeColonyExtent();

        if (I_max == -1) {
            cerr << "\t>>No cells to start 'patient zero' infection! Exiting..<<" << endl;
            f_log << ">>No cells to start 'patient zero' infection! Exiting..<<" << endl;
            exit = true;
        } else {
            cells[I_max][4]++;
            M_tot++;
            numI++;
            numB--;
        }

        // Make sure no more "patient zeros" are spawned
        singleInfectedCell = false;

    }
    if ((uniformPhageInvasion) and (P_0 > 0)) {

        if (not (reflectingBoundary or absorbingBoundary)) {
            cerr << "\t>>Can only do uniform phage invasion when using the box boundary schemes! Exiting..<<" << endl;
            f_log << ">>Can only do uniform phage invasion when using the box boundary schemes! Exiting..<<" << endl;
            exit = true;
        }

        // Update the colony extent
        ComputeColonyExtent();

        // Spawn phages uniformly within the space
        int M = (int)(P_0*pow(L_box,3));
        phages.resize(M);
        for (int m = 0; m < M; m++) {

            double x = 0;
            double y = 0;
            double z = 0;

            while (sqrt(pow(x,2)+pow(y,2)+pow(z,2)) <= r_max) {
                x = (rand(rng)-0.5)*L_box;
                y = (rand(rng)-0.5)*L_box;
                z = (rand(rng)-0.5)*L_box;
            }

            // Spawn the new phage
            phages[m] = vector<double> {x, y, z};

        }

        // Make sure no more phages are artificially spawned
        uniformPhageInvasion = false;

    }
    if (manyInfectedCells) {

        if (not (reflectingBoundary or absorbingBoundary)) {
            cerr << "\t>>Can only do uniform phage invasion when using the box boundary schemes! Exiting..<<" << endl;
            f_log << ">>Can only do uniform phage invasion when using the box boundary schemes! Exiting..<<" << endl;
            exit = true;
        }

        // Update the colony extent
        ComputeColonyExtent();

        // Infect the surface at a depth of dR
        for (int n = 0; n < cells.size(); n++) {
            double r = sqrt(pow(cells[n][0]-center[0],2) + pow(cells[n][1]-center[1],2) + pow(cells[n][2]-center[2],2));
            if (r_max-r < dR) {

                // Choose infection pathway
                if (rand(rng) <= alpha) {  // Lysogenic
                    cells[n][4]--;
                    numL++;
                    numB--;
                } else {                // Lytic
                    cells[n][4] = RandI(9)+1;
                    numI++;
                    numB--;
                }
            }
        }

        // Make sure no more cells are artificially infected
        manyInfectedCells = false;

        // Spawn phages uniformly within the space
        int M = (int)(P_0*pow(L_box,3));
        phages.resize(M);
        for (int m = 0; m < M; m++) {

            double x = 0;
            double y = 0;
            double z = 0;

            while (sqrt(pow(x,2)+pow(y,2)+pow(z,2)) <= r_max) {
                x = (rand(rng)-0.5)*L_box;
                y = (rand(rng)-0.5)*L_box;
                z = (rand(rng)-0.5)*L_box;
            }

            // Spawn the new phage
            phages[m] = vector<double> {x, y, z};
        }
    }
    if (forcedThickness > 0) {

        ComputeColonyExtent();

        int N = cells.size();
        for (int n = 0; n < N; n++) {
            for (int j = 0; j < 3; j++) {
                center[j] += 1/N * cells[n][j];
            }
        }


        // Keep infecting the outermost non-infected cell until number is reached.
        for (int i = 0; i < forcedThickness; i++) {
            r_max = 0;
            I_max = -1;
            for (int n = 0; n < N; n++) {

                if (cells[n][4] > 0) continue;

                double r = sqrt(pow(cells[n][0]-center[0],2) + pow(cells[n][1]-center[1],2) + pow(cells[n][2]-center[2],2));
                if (r >= r_max) {
                    r_max = r;
                    I_max = n;
                }
            }

            cells[I_max][4] = 1;
            numB--;
            numI++;
        }

        forcedThickness = 0;
    }
}

// Replaces the most central bacteria with abstract biomass
void Simulation::CreateCentralBiomass() {

    if (Time*dt < T_i) {return;}
    if (cells.size() == 0) {return;}
    if (numB < 250) {return;}

    // Update the colony extent
    ComputeColonyExtent();

    // Detect which grid point center belongs to
    int i = round( center[0] / dGrid ) + res;
    int j = round( center[1] / dGrid ) + res;
    int k = round( center[2] / dGrid ) + res;

    // Determine the most central bacteria
    double r_min = HUGE_VAL;
    int I = -1;
    vector<int> NN = NearestNeighbours(i,j,k);
    for (int n = 0; n < NN.size(); n++) {

        if (cells[n][4]==0) {
            double r = sqrt(pow(cells[NN[n]][0]-center[0],2) + pow(cells[NN[n]][1]-center[1],2) + pow(cells[NN[n]][2]-center[2],2));

            if (r < r_min) {
                r_min = r;
                I = NN[n];
            }
        }
    }

    if (I == -1) { return; }

    // Compute the maximal distance from central bacteria (I) to any infected bacteria
    r_min = HUGE_VAL;
    for (int n = 0; n < cells.size(); n++) {

        if (cells[n][4]<= 0) continue;

        double r = sqrt(pow(cells[n][0]-cells[I][0],2) + pow(cells[n][1]-cells[I][1],2) + pow(cells[n][2]-cells[I][2],2));

        if (r < r_min) {
            r_min = r;
        }
    }

    // If no cells are yet infected, set r_min to zero to prevet biomass creation
    if (r_min == HUGE_VAL) {r_min = 0.0;}

    // Check if this distance is large enough that the biomass should be created
    if (r_min - cells[I][3]/pow(0.64,1/3) > 1.1*dShell) {

        if (not cells[I][4] == 0) return;

        // Replace the most central bacteria with the central mass
        biomass.resize(8);
        for (int i = 0; i < 8; i++) {
            biomass[i] = cells[I][i];
        }

        // Remove cell from cell list
        cells.erase(cells.begin() + I);

        // Update the CN grid
        UpdateNearestNeighbourGrid();
    }
}


// Returns the gradient of the potential
double Simulation::PotentialGradient(double d, double r1, double r2) {

    // Potential used is harmonic
    // return k*(r1+r2-d)/pow(r1+r2,2);
    return k*(r1+r2-d);
}


// Returns the probabilty for a phage to enter computational regio
double Simulation::FirstPassageTimeDistribution() {

    // Compute the relevant distance to plug into the distribution
    double d = 0.0;
    d = h_agar - h_cell - ( res + 0.5 ) * dGrid;


    // Compute the probability of particle passing the distance at current time
    double p = 0.0;

    for (int k = 1; k < 500; k += 4) {
        p += D_P * (   k  *M_PIl)/pow(d,2) * exp( - pow( (   k  *M_PIl )/(2*d),2) * D_P * Time * dt);
        p -= D_P * ( (k+2)*M_PIl)/pow(d,2) * exp( - pow( ( (k+2)*M_PIl )/(2*d),2) * D_P * Time * dt);
    }

    // Ensure result is positive
    if (p<0) {
        p = 0.0;
    }

    return p;
}


// Computes the center and extent of the colony
void Simulation::ComputeColonyExtent() {

    assert(center.size() == 4);
    if ((int)round(center[3]) == Time) {
        return;
    } else {
        center[0] = 0;
        center[1] = 0;
        center[2] = 0;
        center[3] = Time;
    }

    int N = cells.size();
    for (int n = 0; n < N; n++) {
        for (int j = 0; j < 3; j++) {
            center[j] += 1/N * cells[n][j];
        }
    }
    // Maximal distance from the center
    r_max = 0;
    I_max = -1;
    for (int n = 0; n < N; n++) {
        double r = sqrt(pow(cells[n][0]-center[0],2) + pow(cells[n][1]-center[1],2) + pow(cells[n][2]-center[2],2));
        if (r >= r_max) {
            r_max = r;
            I_max = n;
        }
    }

    // Add safety buffer
    r_max += R;
}


// Updates the CN grid if cells have been removed
void Simulation::UpdateNearestNeighbourGrid() {

    // Update the CN grid
    time_t timer;
    time_t subtimer;
    time(&timer);
    for (int x = 0; x < 2*res+1; x++ ) {
        for (int y = 0; y < 2*res+1; y++ ) {
            for (int z = 0; z < 2*res+1; z++) {
                if (not CN[x][y][z].empty()) {
                    CN[x][y][z].clear();
                }
            }
        }
    }
    benchMark[16] += difftime(time(NULL),timer);
    benchMark[11] += difftime(time(NULL),timer);

    for (int n = 0; n < cells.size(); n++) {

        // Add cell to clever neighbour grid
        time(&subtimer);
        CN[cells[n][5]][cells[n][6]][cells[n][7]].push_back(n);

        benchMark[17] += difftime(time(NULL),subtimer);
        benchMark[11] += difftime(time(NULL),subtimer);
    }

    benchMark[5] += difftime(time(NULL),timer);
}


// Returns all cells in the vicinity of gridpoint (i,j,k)
std::vector<int> Simulation::NearestNeighbours(int i,int j,int k) {

    // Define return vector
    std::vector<int> NN;
    NN.reserve((int)round(pow(3*dGrid,3)));

    // Find neighbouring gridpoints
    for (int ii = i-1; ii <= i+1; ii++ ) {                      // Loop over dim 1
        for (int jj = j-1; jj <= j+1; jj++ ) {                  // Loop over dim 2
            for (int kk = k-1; kk <= k+1; kk++ ) {              // Loop over dim 3
                if ( (ii >= 0) and (ii < 2*res+1) ) {           // Boundary check, dim 1
                    if ( (jj >= 0) and (jj < 2*res+1) ) {       // Boundary check, dim 2
                        if ( (kk >= 0) and (kk < 2*res+1) ) {   // Boundary check, dim 3

                            // Append the CN list to NN
                            NN.insert(NN.end(), CN[ii][jj][kk].begin(), CN[ii][jj][kk].end());

                        }
                    }
                }
            }
        }
    }

    return NN;
}


// Settings /////////////////////////////////////////////////////////////////////////////
// Set the size of the time-step
void Simulation::TimeStep(double dT) {
    this->dt = dT;
    this->dT = dT;
}


// Set the maximum "step" size
void Simulation::TimeStepSkip(int step) {
    maxStep = step;
    posibleStep = step;
}

// Set the maximum diffusion length (in sigma)
void Simulation::MaxStep(double sigma){this->sigma=sigma;}


// Set the type of boundary conditios (1 = simple, 2 = experiment)
void Simulation::BoundaryType(int type) {

    if (type == 1) {
        reflectingBoundary = false;
        absorbingBoundary  = true;
        experimentBoundary = false;
    } else if (type == 2) {
        reflectingBoundary = true;
        absorbingBoundary  = false;
        experimentBoundary = false;
    } else if (type == 3) {
        reflectingBoundary = false;
        absorbingBoundary  = false;
        experimentBoundary = true;
    } else {
        cerr << "\t>>Boundary type not recognized! Exiting..<<" << endl;
        f_log << ">>Boundary type not recognized! Exiting..<<" << endl;
        exit = true;
    }
}


// Set the size of the clever neighbour grid
void Simulation::GridSize(double dGrid){this->dGrid = dGrid;}


// Set the size of the shell margin when using central biomass
void Simulation::ShellSize(double dShell){this->dShell = dShell;}


// Sets the maximum growthrate
void Simulation::MaxGrowthRate(double g_max) {this->g_max = g_max;}


// Sets the maximum colonysize
void Simulation::MaxColonySize(int N_max) {this->N_max = N_max;}


// Sets the time when the phages should start infecting
void Simulation::PhageInvasionStartTime(double T_i) {this->T_i = T_i;}


// Set the type of phage invasion (1 = single, 2 = planar, 3 = uniform, 4 = shell)
void Simulation::PhageInvasionType(int type) {

    invasionType = type;

    if (type == 1) {
        singleInfectedCell   = true;
        planarPhageInvasion  = false;
        uniformPhageInvasion = false;
        manyInfectedCells    = false;
        forcedThickness      = 0;
    } else if (type == 2) {
        singleInfectedCell   = false;
        planarPhageInvasion  = true;
        uniformPhageInvasion = false;
        manyInfectedCells    = false;
        forcedThickness      = 0;
    } else if (type == 3) {
        singleInfectedCell   = false;
        planarPhageInvasion  = false;
        uniformPhageInvasion = true;
        manyInfectedCells    = false;
        forcedThickness      = 0;
    } else if (type == 4) {
        singleInfectedCell   = false;
        planarPhageInvasion  = false;
        uniformPhageInvasion = false;
        manyInfectedCells    = true;
        forcedThickness      = 0;
    } else if (type == 5) {
        singleInfectedCell   = true;
        planarPhageInvasion  = false;
        uniformPhageInvasion = true;
        manyInfectedCells    = false;
        forcedThickness      = 0;
    }
}
// Set the type of phage invasion (1 = single, 2 = planar, 3 = uniform, 4 = shell)
void Simulation::PhageInvasionType(int type, int thickness) {

    invasionType = type;

    if (type == 1) {
        singleInfectedCell   = true;
        planarPhageInvasion  = false;
        uniformPhageInvasion = false;
        manyInfectedCells    = false;
        forcedThickness      = thickness;
    } else if (type == 2) {
        singleInfectedCell   = false;
        planarPhageInvasion  = true;
        uniformPhageInvasion = false;
        manyInfectedCells    = false;
        forcedThickness      = thickness;
    } else if (type == 3) {
        singleInfectedCell   = false;
        planarPhageInvasion  = false;
        uniformPhageInvasion = true;
        manyInfectedCells    = false;
        forcedThickness      = thickness;
    } else if (type == 4) {
        singleInfectedCell   = false;
        planarPhageInvasion  = false;
        uniformPhageInvasion = false;
        manyInfectedCells    = true;
        forcedThickness      = thickness;
    } else if (type == 5) {
        singleInfectedCell   = true;
        planarPhageInvasion  = false;
        uniformPhageInvasion = true;
        manyInfectedCells    = false;
        forcedThickness      = thickness;
    } else if (type == 6) {
        singleInfectedCell   = false;
        planarPhageInvasion  = false;
        uniformPhageInvasion = false;
        manyInfectedCells    = false;
        forcedThickness      = thickness;
    }

}

// Sets initial density of the phages (1/µm^2)
void Simulation::PhageInitialDensity(double P_0) {
    if (P_0 >= 0) {
        this->P_0 = P_0;
    }
}


// Sets the phage parameters to De Paepe values for "phageType"
void Simulation::PhageType(std::string& phageType) {

    // Load P1 parameters
    if (!phageType.compare("P1vir")) {
        gamma = 1/dT;                               //          Probability to infect cell
        alpha = 0.00;                               //          Probability for phage to go lysogenic
        beta  = 400;                                //          Multiplication factor phage
        delta = 0.0032;                             // [1/hour] Rate of phage decay
        r     = 10.0/(60.0/60.0);                   //          Constant used in the time-delay mechanism
        r    *= (1+K/(n_0*pow(dGrid,3)))/g_max;
        D_P   = 4000;                               // [µm^2/hour] Diffusion constant for the phage
        this->phageType = phageType;
    }

    // Load lambda parameters
    if (!phageType.compare("lambdavir")) {
        gamma = 1/dT;                               //          Probability to infect cell
        alpha = 0.00;                               //          Probability for phage to go lysogenic
        beta  = 115;                                //          Multiplication factor phage
        delta = 0.0030;                             // [1/hour] Rate of phage decay
        r     = 10.0/(42.0/60.0);                   //          Constant used in the time-delay mechanism
        r    *= (1+K/(n_0*pow(dGrid,3)))/g_max;
        D_P   = 40000;                              // [µm^2/hour] Diffusion constant for the phage
        this->phageType = phageType;
    }

    // Load P1 parameters
    if (!phageType.compare("T4")) {
        gamma = 1/dT;                               //          Probability to infect cell
        alpha = 0.00;                               //          Probability for phage to go lysogenic
        beta  = 150;                                //          Multiplication factor phage
        delta = 0.00283;                            // [1/hour] Rate of phage decay
        r     = 10.0/(23.0/60.0);                   //          Constant used in the time-delay mechanism
        r    *= (1+K/(n_0*pow(dGrid,3)))/g_max;
        D_P   = 13000;                              // [µm^2/hour] Diffusion constant for the phage
        this->phageType = phageType;
    }
}


// Sets the diffusion constant of the phages
void Simulation::PhageDiffusionConstant(double D_P) {
    this->D_P = D_P;
    phageType = "user";
}


// Sets rate of the infection increaasing in stage
void Simulation::PhageInfectionRate(double r) {
    this->r = 10*r;
    phageType = "user";
}


// Set the decay rate of the phages
void Simulation::PhageDecayRate(double delta) {
    this->delta = delta;
    phageType = "user";
}


// Set the lysogen frequency of the phage
void Simulation::PhageLysogenFrequency(double alpha) {
    this->alpha = alpha;
    phageType = "user";
}


// Set the size of the bursts
void Simulation::PhageBurstSize(int beta) {
    this->beta = beta;
    phageType = "user";
}


// Changes the adsorption parameter gamma
void Simulation::PhageAdsorptionParameter(double gamma) {
    this->gamma = gamma;
    phageType = "user";
}


// Changes the mutation/resistance rate of the cells
void Simulation::CellMutationRate(double epsilon) {this->epsilon = epsilon;}


// Sets initial density of the bacteria (1/µm^2)
void Simulation::CellInitialDensity(double B_0) {
    if (B_0 >= 0) {
        L_box = 1/sqrt(B_0);
        res = floor(L_box/(2*dGrid)-1/2);
    }
}


// Helping functions ////////////////////////////////////////////////////////////////////
// Returns random integter between 0 and n
int Simulation::RandI(int n) {

    // Set limit on distribution
    uniform_int_distribution <int> distr(0, n);

    return distr(rng);
}


// Returns random double between 0 and n
double Simulation::Rand(double n) {

    // Set limit on distribution
    uniform_real_distribution <double> distr(0, n);

    return distr(rng);
}


// Returns random normal dist. number with mean m and variance s^2
double Simulation::RandN(double m, double s) {

    // Set limit on distribution
    normal_distribution <double> distr(m, s);

    return distr(rng);
}


// Returns the current number of uninfected cells
int Simulation::NumberOfUninfectedCells() {
    return numB;
}


// Returns the current number of phages
int Simulation::NumberOfPhages() {
    return phages.size();
}


// Returns the total volume of all the cells
double Simulation::Biomass() {
    double V = 0.0;
    for (int n = 0; n < cells.size(); n++) {
        V += 4*M_PIl/3*pow(cells[n][3],3);
    }
    return V;
}

// Returns the current gamma value
double Simulation::GetGamma() {
    return gamma;
}

// Returns the current r_max
double Simulation::GetColonyExtent() {
    ComputeColonyExtent();
    return r_max;
}

// Sets the seed of the random number generator
void Simulation::SetRngSeed(int n) {
    rngSeed = n;
    rng.seed( rngSeed );
}


// Debug function (prints "debug" and int)
void Simulation::deb(int n) {

    std::stringstream stream;
    stream << "Debugging: " << n << endl;
    cout << stream.str();
}


// Write a log.txt file
void Simulation::WriteLog(double T) {
    if ((not f_log.is_open()) and (not exit)) {

        // Open the file stream and write the command
        f_log.open(path + "/log.txt", fstream::out);

        // Physical parameters
        f_log << "N_0 = "      << N_0      << endl;        // Number of starting cells
        f_log << "P_0 = "      << P_0      << endl;        // Density of invading phages
        f_log << "N_max = "    << N_max    << endl;        // Maximal allowed cells
        f_log << "M_max = "    << M_max    << endl;        // Maximum allowed phages
        f_log << "dR = "       << dR       << endl;        // Depth of the infected layer (for spawning)
        f_log << "dShell = "   << dShell   << endl;        // Thickness of agent layer (for biomass)
        f_log << "n_0 = "      << n_0      << endl;        // Initial nutrient density
        f_log << "g_max = "    << g_max    << endl;        // Maximal growth rate
        f_log << "K = "        << K        << endl;        // Michaels-Menten constant for Monod growth
        f_log << "r = "        << r        << endl;        // Constant used in the time-delay mechanism
        f_log << "T_i = "      << T_i      << endl;        // Phage infection time
        f_log << "L_box = "    << L_box    << endl;        // Length of boundary condition box
        f_log << "k = "        << k        << endl;        // Strength of the repulsion potential
        f_log << "R = "        << R        << endl;        // Critical radius
        f_log << "D_B = "      << D_B      << endl;        // Bacteria difusion
        f_log << "D_P = "      << D_P      << endl;        // Phage diffusion
        f_log << "D_n = "      << D_n      << endl;        // Nutrient diffusion
        f_log << "phageType = "<< phageType<< endl;        // Name of the phage being used
        f_log << "alpha = "    << alpha    << endl;        // Probability for phage to go lysogenic
        f_log << "beta = "     << beta     << endl;        // Multiplication factor phage
        f_log << "gamma = "    << gamma    << endl;        // Probability to infect cell
        f_log << "delta = "    << delta    << endl;        // Rate of phage decay
        f_log << "epsilon = "  << epsilon  << endl;        // Rate of resistance mutations
        f_log << "eta = "      << eta      << endl;        // Amount of division noise
        f_log << "nu = "       << nu       << endl;        // Amount of displacement noise

        f_log << "h_agar = "   << h_agar   << endl;        // Total height of the soft agar layers
        f_log << "h_cell = "   << h_cell   << endl;        // Distance of the cell colony to the hard agar

        // Non physical parameters
        f_log << "res = "      << res      << endl;        // Resolution of grid
        f_log << "dGrid = "    << dGrid    << endl;        // Spacing of the grid
        f_log << "dt = "       << dt       << endl;        // Default time step
        f_log << "dT_c = "     << dT_c     << endl;        // Time step size for cell updates

        f_log << "maxStep = "  << maxStep  << endl;        // Time step step size
        f_log << "sigma = "    << sigma    << endl;        // Maximum number of standard deviations to diffuse
        f_log << "rngSeed = "  << rngSeed  << endl;        // Random Number Generator seed

        f_log << "singleInfectedCell = "   << singleInfectedCell   << endl;    //
        f_log << "planarPhageInvasion = "  << planarPhageInvasion  << endl;    //
        f_log << "uniformPhageInvasion = " << uniformPhageInvasion << endl;    // Phage Invasion type
        f_log << "manyInfectedCells = "    << manyInfectedCells    << endl;    //
        f_log << "forcedThickness = "      << forcedThickness      << endl;    //

        f_log << "reflectingBoundary = "   << reflectingBoundary   << endl;    //
        f_log << "absorbingBoundary = "    << absorbingBoundary    << endl;    // Boundary type
        f_log << "experimentBoundary = "   << experimentBoundary   << endl;    //
    }
}


// Store the state of the simulation
void Simulation::SaveState() {

    // Open export stream
    std::ofstream f_out;

    // Save the state of the simulation to a file
    f_out.open(path + "/state/rng.txt",fstream::out);
    f_out << rng;
    f_out.close();

    // Flatten cells matrix into cells_flat vector
    std::vector<double> cells_flat;
    if (not biomass.empty()) {
      for (auto el : biomass) {
        cells_flat.push_back(el);
      }
    }
    for (auto vec : cells) {
      for (auto el : vec) {
        cells_flat.push_back(el);
      }
    }
    // Convert flattened vector into arma matrix
    mat C(cells_flat);
    C = reshape(C,8,C.n_elem/8).t();
    C.save(path + "/state/cells.txt");

    // Flatten phage matrix into phage_flat vector
    std::vector<double> phage_flat;
    for (auto vec : phages) {
      for (auto el : vec) {
        phage_flat.push_back(el);
      }
    }
    // Convert flattened vector into arma matrix
    mat P(phage_flat);
    P = reshape(P,8,P.n_elem/8).t();
    P.save(path + "/state/phages.txt");

    // Export nutrient level
    if (nutrientField) {
        nutrient_grid.save(path + "/state/nutrient.txt");
    } else {
        std::ofstream f_out(path + "/state/nutrient.txt",fstream::out);
        f_out.precision(16);
        f_out << fixed << nutrient;
        f_out.close();
    }

    // Export the simulation settngs
    f_out.open(path + "/state/settings.txt",fstream::out);
    f_out.precision(16);

    // Simulation parameters
    f_out << "N_0="                     << N_0                      << endl;    // Number of starting cells
    f_out << "P_0="     << fixed        << P_0                      << endl;    // Density of invading phages

    f_out << "L_box="   << fixed        << L_box                    << endl;    // Length of boundary condition box
    f_out << "dR="      << fixed        << dR                       << endl;    // Depth of the infected layer (for spawning)
    f_out << "dShell="  << fixed        << dShell                   << endl;    // Thickness of agent layer (for biomass)

    f_out << "h_agar="  << fixed        << h_agar                   << endl;    // Total height of the soft agar layers
    f_out << "h_cell="  << fixed        << h_cell                   << endl;    // Distance of the cell colony to the hard agar

    // Simulation tallies
    f_out << "M_tot="                   << M_tot                    << endl;    // The accumulated number of phages
    f_out << "numB="                    << numB                     << endl;    // Current tally of suceptible cells
    f_out << "numL="                    << numL                     << endl;    // Current tally of lysogenic cells
    f_out << "numI="                    << numI                     << endl;    // Current tally of infected cells
    f_out << "numD="                    << numD                     << endl;    // Current tally of dead cells
    f_out << "numR="                    << numR                     << endl;    // Current tally of resistant cells

    // Physical parameters
    f_out << "n_0="     << fixed        << n_0                      << endl;    // Initial nutrient density
    f_out << "g_max="   << fixed        << g_max                    << endl;    // Maximal growth rate
    f_out << "K="       << fixed        << K/(n_0*pow(dGrid,3))     << endl;    // Michaels-Menten constant for Monod growth
    f_out << "r="       << fixed        << r                        << endl;    // Constant used in the time-delay mechanism
    f_out << "T_i="     << fixed        << T_i                      << endl;    // Phage infection time
    f_out << "k="       << fixed        << k                        << endl;    // Strength of the repulsion potential
    f_out << "R="       << fixed        << R                        << endl;    // Division radius

    f_out << "D_B="     << fixed        << D_B                      << endl;    // Bacteria difusion
    f_out << "D_P="     << fixed        << D_P                      << endl;    // Phage diffusion
    f_out << "D_n="     << fixed        << D_n                      << endl;    // Nutrient diffusion

    f_out << "alpha="   << fixed        << alpha                    << endl;    // Probability for phage to go lysogenic
    f_out << "beta="                    << beta                     << endl;    // Multiplication factor phage
    f_out << "gamma="   << fixed        << gamma                    << endl;    // Probability to infect cell
    f_out << "delta="   << fixed        << delta                    << endl;    // Rate of phage decay
    f_out << "epsilon= "<< fixed        << epsilon                  << endl;    // Rate of resistance mutations

    f_out << "eta="     << fixed        << eta                      << endl;    // Amount of division noise
    f_out << "nu="      << fixed        << nu                       << endl;    // Amount of displacement noise

    f_out << "phageType="               << phageType                << endl;    // Contains the type of phage chosen (sets parameters according to De Paepe)
    f_out << "rngSeed=" << fixed        << rngSeed                  << endl;    // Initial seed for the random number generator

    // Non physical parameters
    f_out << "Time="                    << Time                     << endl;    // Integer to keep track of time
    f_out << "dt="      << fixed        << dt                       << endl;    // Default Time step size
    f_out << "dT_c="    << fixed        << dT_c                     << endl;    // Time step size for cell updates
    f_out << "dT="      << fixed        << dT                       << endl;    // Time step size for everything else

    f_out << "maxStep="                 << maxStep                  << endl;    // Maximum time-steps to "step" when finegrained bool is false
    f_out << "sigma="                   << sigma                    << endl;    // Maximum number of standard deviations to move in diffusions term
    f_out << "posibleStep="             << posibleStep              << endl;    // Largest possible "step" at current time

    f_out << "N_max="                   << N_max                    << endl;    // Maximal allowed cells
    f_out << "M_max="                   << M_max                    << endl;    // Maximum allowed phages

    f_out << "res="                     << res                      << endl;    // Resolution of grid
    f_out << "dGrid="   << fixed        << dGrid/R                  << endl;    // Spacing of the grid

    f_out << "p="       << fixed        << p                        << endl;    // Helper variable to control phage spawning

    f_out << "debug="                   << debug                    << endl;    // Set the debug level
    f_out << "invasionType="            << invasionType             << endl;    // Number to remember the chose invasion type.

    f_out << "path="                    << path                     << endl;    // Sets the path to store in

    // Booleans
    f_out << boolalpha;
    f_out << "exit="                    << exit                     << endl;    // Exit bool; terminates execution
    f_out << "debugBool="               << debugBool                << endl;

    f_out << "nutrientField="           << nutrientField            << endl;    // Boolean to toggle between nutrient field and just nutrients

    f_out << "exportAny="               << exportAny                << endl;    //
    f_out << "exportCellData="          << exportCellData           << endl;    //
    f_out << "exportColonySize="        << exportColonySize         << endl;    //
    f_out << "exportPhageData="         << exportPhageData          << endl;    // Booleans to control the export output
    f_out << "exportNutrient="          << exportNutrient           << endl;    //

    f_out << "singleInfectedCell="      << singleInfectedCell       << endl;    // Boolean to enable phages via a single infected cell
    f_out << "planarPhageInvasion="     << planarPhageInvasion      << endl;    // Boolean to enable phages invading in a plane
    f_out << "uniformPhageInvasion="    << uniformPhageInvasion     << endl;    // Boolean to enable phages invading uniformly from entire space
    f_out << "manyInfectedCells="       << manyInfectedCells        << endl;    // Boolean to enable phages invading from an infected surface and uniformly from entire space
    f_out << "forcedThickness="         << forcedThickness          << endl;    // "Boolean" to force an infected surface of a given number

    f_out << "reflectingBoundary="      << reflectingBoundary       << endl;    // Boolean to enable the simple box boundary condition with reflective sides
    f_out << "absorbingBoundary="       << absorbingBoundary        << endl;    // Boolean to enable the simple box boundary condition with absorbing sides
    f_out << "experimentBoundary="      << experimentBoundary       << endl;    // Boolean to enable the advanced petri dish boundary conditions

    f_out.close();
}


// Load a previous state
void Simulation::LoadState(string& loadPath) {

    loadPath = "data/" + loadPath;

    // Verify that all files exist
    struct stat info;
    string testPath = loadPath;
    if (not (stat(testPath.c_str(), &info) == 0     && S_ISDIR(info.st_mode))) { exit = true; }
    testPath = loadPath + "/state";
    if (not (stat(testPath.c_str(), &info) == 0     && S_ISDIR(info.st_mode))) { exit = true; }
    testPath = loadPath + "/state/settings.txt";
    if (not (stat(testPath.c_str(), &info) == 0     && S_ISREG(info.st_mode))) { exit = true; }
    testPath = loadPath + "/state/rng.txt";
    if (not (stat(testPath.c_str(), &info) == 0     && S_ISREG(info.st_mode))) { exit = true; }
    testPath = loadPath + "/state/cells.txt";
    if (not (stat(testPath.c_str(), &info) == 0     && S_ISREG(info.st_mode))) { exit = true; }
    testPath = loadPath + "/state/phages.txt";
    if (not (stat(testPath.c_str(), &info) == 0     && S_ISREG(info.st_mode))) { exit = true; }
    testPath = loadPath + "/state/nutrient.txt";
    if (not (stat(testPath.c_str(), &info) == 0     && S_ISREG(info.st_mode))) { exit = true; }
    if (exit) {
        cout << "\t>>Could not verify loadPath! Exiting..<<" << endl;
        return;
    }

    // Load settings
    ifstream is_file(loadPath + "/state/settings.txt");

    string line;
    while( getline(is_file, line) ) {

        istringstream is_line(line);
        string key;

        if( getline(is_line, key, '=') ) {
            string value;

            if( getline(is_line, value) ) {

                // Simulation parameters
                     if (key.compare("N_0")                    == 0) { istringstream(value) >> N_0; }
                else if (key.compare("P_0")                    == 0) { istringstream(value) >> P_0; }
                else if (key.compare("L_box")                  == 0) { istringstream(value) >> L_box; }
                else if (key.compare("dR")                     == 0) { istringstream(value) >> dR; }
                else if (key.compare("dShell")                 == 0) { istringstream(value) >> dShell; }
                else if (key.compare("h_agar")                 == 0) { istringstream(value) >> h_agar; }
                else if (key.compare("h_cell")                 == 0) { istringstream(value) >> h_cell; }
                // Simulation tallies
                else if (key.compare("M_tot")                  == 0) { istringstream(value) >> M_tot; }
                else if (key.compare("numB")                   == 0) { istringstream(value) >> numB; }
                else if (key.compare("numL")                   == 0) { istringstream(value) >> numL; }
                else if (key.compare("numI")                   == 0) { istringstream(value) >> numI; }
                else if (key.compare("numD")                   == 0) { istringstream(value) >> numD; }
                else if (key.compare("numR")                   == 0) { istringstream(value) >> numR; }
                // Physical parameters
                else if (key.compare("n_0")                    == 0) { istringstream(value) >> n_0; }
                else if (key.compare("g_max")                  == 0) { istringstream(value) >> g_max; }
                else if (key.compare("K")                      == 0) { istringstream(value) >> K; }
                else if (key.compare("r")                      == 0) { istringstream(value) >> r; }
                else if (key.compare("T_i")                    == 0) { istringstream(value) >> T_i; }
                else if (key.compare("k")                      == 0) { istringstream(value) >> k; }
                else if (key.compare("R")                      == 0) { istringstream(value) >> R; }
                else if (key.compare("D_B")                    == 0) { istringstream(value) >> D_B; }
                else if (key.compare("D_P")                    == 0) { istringstream(value) >> D_P; }
                else if (key.compare("D_n")                    == 0) { istringstream(value) >> D_n; }
                else if (key.compare("alpha")                  == 0) { istringstream(value) >> alpha; }
                else if (key.compare("beta")                   == 0) { istringstream(value) >> beta; }
                else if (key.compare("gamma")                  == 0) { istringstream(value) >> gamma; }
                else if (key.compare("delta")                  == 0) { istringstream(value) >> delta; }
                else if (key.compare("epsilon")                == 0) { istringstream(value) >> epsilon; }
                else if (key.compare("eta")                    == 0) { istringstream(value) >> eta; }
                else if (key.compare("nu")                     == 0) { istringstream(value) >> nu; }
                else if (key.compare("phageType")              == 0) { istringstream(value) >> phageType; }
                else if (key.compare("rngSeed")                == 0) { istringstream(value) >> rngSeed; }
                // Non physical parameters
                else if (key.compare("Time")                   == 0) { istringstream(value) >> Time; }
                else if (key.compare("dt")                     == 0) { istringstream(value) >> dt; }
                else if (key.compare("dT_c")                   == 0) { istringstream(value) >> dT_c; }
                else if (key.compare("dT")                     == 0) { istringstream(value) >> dT; }
                else if (key.compare("maxStep")                == 0) { istringstream(value) >> maxStep; }
                else if (key.compare("sigma")                  == 0) { istringstream(value) >> sigma; }
                else if (key.compare("posibleStep")            == 0) { istringstream(value) >> posibleStep; }
                else if (key.compare("N_max")                  == 0) { istringstream(value) >> N_max; }
                else if (key.compare("M_max")                  == 0) { istringstream(value) >> M_max; }
                else if (key.compare("res")                    == 0) { istringstream(value) >> res; }
                else if (key.compare("dGrid")                  == 0) { istringstream(value) >> dGrid; }
                else if (key.compare("p")                      == 0) { istringstream(value) >> p; }
                else if (key.compare("debug")                  == 0) { istringstream(value) >> debug; }
                else if (key.compare("invasionType")           == 0) { istringstream(value) >> invasionType; }
                else if (key.compare("path")                   == 0) { istringstream(value) >> path; }
                // Booleans
                else if (key.compare("exit")                   == 0) { istringstream(value) >> boolalpha >> exit; }
                else if (key.compare("debugBool")              == 0) { istringstream(value) >> boolalpha >> debugBool; }
                else if (key.compare("nutrientField")          == 0) { istringstream(value) >> boolalpha >> nutrientField; }
                else if (key.compare("exportAny")              == 0) { istringstream(value) >> boolalpha >> exportAny; }
                else if (key.compare("exportCellData")         == 0) { istringstream(value) >> boolalpha >> exportCellData; }
                else if (key.compare("exportColonySize")       == 0) { istringstream(value) >> boolalpha >> exportColonySize; }
                else if (key.compare("exportPhageData")        == 0) { istringstream(value) >> boolalpha >> exportPhageData; }
                else if (key.compare("exportNutrient")         == 0) { istringstream(value) >> boolalpha >> exportNutrient; }
                else if (key.compare("singleInfectedCell")     == 0) { istringstream(value) >> boolalpha >> singleInfectedCell; }
                else if (key.compare("planarPhageInvasion")    == 0) { istringstream(value) >> boolalpha >> planarPhageInvasion; }
                else if (key.compare("uniformPhageInvasion")   == 0) { istringstream(value) >> boolalpha >> uniformPhageInvasion; }
                else if (key.compare("manyInfectedCells")      == 0) { istringstream(value) >> boolalpha >> manyInfectedCells; }
                else if (key.compare("forcedThickness")        == 0) { istringstream(value) >> forcedThickness; }
                else if (key.compare("reflectingBoundary")     == 0) { istringstream(value) >> boolalpha >> reflectingBoundary; }
                else if (key.compare("absorbingBoundary")      == 0) { istringstream(value) >> boolalpha >> absorbingBoundary; }
                else if (key.compare("experimentBoundary")     == 0) { istringstream(value) >> boolalpha >> experimentBoundary; }
            }
        }
    }

    // Initialize the simulation matrices
    Initialize();

    // Load the random number generator state
    ifstream f_in(loadPath + "/state/rng.txt",fstream::in);
    f_in >> rng;
    f_in.close();

    // Load the cell configuration
    mat C;
    C.load(loadPath + "/state/cells.txt");
    cells.resize(C.n_rows);
    for (size_t i = 0; i < C.n_rows; ++i) {
        cells[i] = conv_to< vector<double> >::from(C.row(i));
    };
    // Check if central biomass is used
    if (cells[0][4] > R ) {
        biomass = cells[0];
        cells.erase(cells.begin());
    }

    // Load the phage configuration
    mat P;
    P.load(loadPath + "/state/phages.txt");
    phages.resize(P.n_rows);
    for (size_t i = 0; i < P.n_rows; ++i) {
        phages[i] = conv_to< vector<double> >::from(P.row(i));
    };

    // Export nutrient level
    if (nutrientField) {
        nutrient_grid.save(path + "/state/nutrient.txt");
    } else {
        f_in.open(path + "/state/nutrient.txt",fstream::in);
        f_in >> nutrient;
        f_in.close();
    }
}


// Set debug level to 0
void Simulation::Quiet() { debug=0; }


// Sets exit to false;
void Simulation::ClearErrors() { exit = false;}


// Set the number of output samples
void Simulation::SetSamples(int nSamp) {this->nSamp = nSamp;};


// File outputs /////////////////////////////////////////////////////////////////////////
// Sets boolean for export function
void Simulation::ExportCellData()         { exportCellData          = true; exportAny = true; };
void Simulation::ExportColonySize()       { exportColonySize        = true; exportAny = true; };
void Simulation::ExportPhageData()        { exportPhageData         = true; exportAny = true; };
void Simulation::ExportNutrient()         { exportNutrient          = true; exportAny = true; };


// Master function to export the data
void Simulation::ExportData(double t) {

    // Export the data
    if ( exportCellData ) {
        f_ExportCellData(t);
    }
    if ( exportColonySize ) {
        f_ExportColonySize(t);
    }
    if ( exportPhageData ) {
        f_ExportPhageData(t);
    }
    if ( exportNutrient ) {
        f_ExportNutrient(t);
    }
}


// Export the position and size of the cells
void Simulation::f_ExportCellData(double t) {

    // Verify the file stream is open
    string fileName = "CellData";
    OpenFileStream(f_cells, fileName);

    // Output central biomass i used
    if (not biomass.empty()) {
        // Output format:   (T    x    y    z   r   s)
        f_cells << fixed    << setprecision(8);
        f_cells << setw(8)  << t << "\t";
        f_cells << setw(12) << biomass[0] << "\t";
        f_cells << setw(12) << biomass[1] << "\t";
        f_cells << setw(12) << biomass[2] << "\t";
        f_cells << setw(12) << biomass[3] << "\t";
        f_cells << fixed    << setprecision(0);
        f_cells << setw(12) << biomass[4] << "\t";
        f_cells << setw(12) << 0          << endl;
    }

    // Loop over cells in simulation
    for (int i = 0; i < cells.size(); i++) {

        // Output format:   (T    x    y    z   r   s)
        f_cells << fixed    << setprecision(8);
        f_cells << setw(8)  << t << "\t";
        f_cells << setw(12) << cells[i][0] << "\t";
        f_cells << setw(12) << cells[i][1] << "\t";
        f_cells << setw(12) << cells[i][2] << "\t";
        f_cells << setw(12) << cells[i][3] << "\t";
        f_cells << fixed    << setprecision(0);
        f_cells << setw(12) << cells[i][4] << "\t";
        f_cells << setw(12) << cells[i][8] << endl;
    }
}


// Export the position and size of the phages
void Simulation::f_ExportPhageData(double t) {

    // Verify the file stream is open
    string fileName = "PhageData";
    OpenFileStream(f_phages, fileName);

    // Loop over phages in simulation
    for (int i = 0; i < phages.size(); i++) {

        // Output format:   T    x    y    z
        f_phages << fixed    << setprecision(8);
        f_phages << setw(8)  << t << "\t";
        f_phages << setw(12) << phages[i][0] << "\t";
        f_phages << setw(12) << phages[i][1] << "\t";
        f_phages << setw(12) << phages[i][2] << endl;
    }
}


// Export the volume of colony and number of cells.
void Simulation::f_ExportColonySize(double t) {

    // Verify the file stream is open
    string fileName = "ColonySize";
    OpenFileStream(f_colonySize, fileName);

    // Writes the time, total volume, number of cells,
    // number of lytic stage cells, number of lysogenic cells

    double V    = 0.0;
    int N       = cells.size();
    int M       = phages.size();

    // Loop over cells
    for (int n = 0; n < N; n++) {
        // Increase total volume
        V += 4.0/3.0*M_PIl*pow(cells[n][3],3);
    }
    if (not biomass.empty()) {
        V  += 4.0/3.0*M_PIl*pow(biomass[3],3);
    }

    f_colonySize << fixed    << setprecision(3);
    f_colonySize << setw(6)  << t       << "\t";
    f_colonySize << setw(8)  << V       << "\t";
    f_colonySize << setw(6)  << N       << "\t";
    f_colonySize << setw(6)  << numB    << "\t";
    f_colonySize << setw(6)  << numL    << "\t";
    f_colonySize << setw(6)  << numI    << "\t";
    f_colonySize << setw(6)  << numD    << "\t";
    f_colonySize << setw(6)  << numR    << "\t";
    f_colonySize << setw(6)  << M       << "\t";
    f_colonySize << setw(8)  << M_tot   << endl;
}


// Export the amount of the nutrient
void Simulation::f_ExportNutrient(double t) {

    // Verify the file stream is open
    string fileName = "Nutrient";
    OpenFileStream(f_n, fileName);

    // Compute the nutrient per µm^3
    double N = 0.0;

    if (nutrientField) {
        N = accu(nutrient_grid)/(pow(2*res+1,3)*pow(dGrid,3));
    } else {
        N = nutrient/pow(dGrid,3);
    }

    f_n << setw(12) << t << "\t";
    f_n << setw(12) << N << endl;
}


// Data Handling ////////////////////////////////////////////////////////////////////////
// Open filstream if not allready opened
void Simulation::OpenFileStream(ofstream& stream, string& fileName) {

    // Check that if file stream is open.
    if ((not stream.is_open()) and (not exit)) {

        // Debug info
        if (debug > 0) {cout << "\t\tSaving data to file: " << path << "/" << fileName << ".txt" << endl;}

        // Check if the output file exists
        bool fileExists = false;
        struct stat info;
        string streamPath;
        streamPath = path + "/" + fileName + ".txt";

        if (stat(streamPath.c_str(), &info) == 0 && S_ISREG(info.st_mode)) {
            fileExists = true;
        }

        // Open the file stream
        if (not fileExists) {
            stream.open(streamPath, fstream::out);
        } else {
            stream.open(streamPath, fstream::app);
        }

        // Check stream is open
        if ((not exit) and (not stream.is_open())) {
            cerr << "\t>>Could not open filestream \"" << streamPath << "\"! Exiting..<<" << endl;
            f_log <<  ">>Could not open filestream \"" << streamPath << "\"! Exiting..<<" << endl;
            exit = true;
        };

        // If appending to existing file, do not rewrite the meta data
        if (fileExists) {
            return;
        }

        // Write meta data to the data file
        stream << "Datatype: "  << fileName << endl;
    }
}


// Generates a save path for datafiles
string Simulation::GeneratePath() {

    firstRun = true;

    // Have you allready made the path?
    if (path.compare(0,4,"data") == 0) {
        firstRun = false;
        return path;
    }

    // Are the filestreams allready open?
    if ( f_cells.is_open() or f_colonySize.is_open() or f_phages.is_open() or f_n.is_open() or f_log.is_open() ) {
        firstRun = false;
        return path;
    }

    // Generate a directory path
    string prefix = "data/";    // Data folder name
    time_t t = time(0);         // Get time now
    struct tm  tstruct;         // And format the date
    tstruct = *localtime(&t);   // as "MNT_DD_YY" for folder name
    char buffer[80];            // Create a buffer to store the date
    strftime(buffer, sizeof(buffer), "%F", &tstruct);  // Store the formated foldername in buffer
    string dateFolder(buffer);

    // Construct path from prefix and date folder
    string path_s = prefix;
    path_s += dateFolder;

    // Check if path exists
    struct stat info;
    if (not (stat(path_s.c_str(), &info) == 0 && S_ISDIR(info.st_mode))) {    // Create path if it does not exist
        mkdir(path_s.c_str(),0700);
    }

    int currentNumerateFolder = 1;
    // Check if user has specified numbered folder
    if (path.empty()) {

        // Loop over folders in date folder, to find current number
        DIR *dir;
        if ((dir = opendir (path_s.c_str())) != NULL) {
            struct dirent *ent;
            while ((ent = readdir (dir)) != NULL) {
                if (ent->d_type == DT_DIR) {
                    // Skip . or ..
                    if (ent->d_name[0] == '.') {continue;}
                    currentNumerateFolder++;        // Increment folder number
                }
            }
            closedir (dir);
        }

    } else {    // User has specified a path

        try {   // Check if path is folder number
            currentNumerateFolder = stoi(path);

        } catch (const std::exception& e) { // Path is full path

            // This path maybe more than one layer deep, so attempt to make it recursively
            int len = path.length();

            string folder = "";
            for (int i = 0; i < len; i++) {

                if (path[i] == '/') {

                    // Append folder to path
                    path_s += "/";
                    path_s += folder;

                    // Make folder
                    if(not (stat(path_s.c_str(), &info) == 0 && S_ISDIR(info.st_mode))) {    // Create path if it does not exist
                        mkdir(path_s.c_str(),0700);
                    }

                    folder = "";    // Reset folder
                } else {
                    folder += path[i];  // Append char to folder name
                }

                if (i == len -1 ) {

                    // Append folder to path
                    path_s += "/";
                    path_s += folder;

                    // Make folder
                    if(not (stat(path_s.c_str(), &info) == 0 && S_ISDIR(info.st_mode))) {    // Create path if it does not exist
                        mkdir(path_s.c_str(),0700);

                        // Create a folder for the state
                        string path_ss = path_s + "/state";
                        mkdir(path_ss.c_str(),0700);

                        firstRun = true;
                    }
                }
            }

            return path_s;
        }
    }

    // Append numerate folder
    path_s += "/";
    path_s += to_string(currentNumerateFolder);

    // Check if path exists
    if(not (stat(path_s.c_str(), &info) == 0 && S_ISDIR(info.st_mode))) {    // Create path if it does not exist
        mkdir(path_s.c_str(),0700);

        // Create a folder for the state
        string path_ss = path_s + "/state";
        mkdir(path_ss.c_str(),0700);

        firstRun = true;
    }

    return path_s;
}


// Sets the folder number (useful when running parralel code)
void Simulation::SetFolderNumber(int number) {path = to_string(number);};

// Sets the folder path (useful when running parralel code)
void Simulation::SetPath(std::string& path) {this->path = path;}

// Get properties ///////////////////////////////////////////////////////////////////////
// Returns the save path
std::string Simulation::GetPath() {
    return path;
}


// Returns the time
int Simulation::GetTime() {
    return Time;
}


// Returns the time-step dT
int Simulation::GetDeltaT() {
    return dT;
}


// Operators ////////////////////////////////////////////////////////////////////////////
// Copy assignment
Simulation& Simulation::operator=(const Simulation& rhs){

    N_0         = rhs.N_0;                                // Store the initial number of cells in the simulation
    P_0         = rhs.P_0;                                // [1/µm^2] The density of invading phages in the simulation initially
    N_max       = rhs.N_max;                              // Maximum number of cells in simulation
    M_max       = rhs.M_max;                              // Maximum number of phages in simulation
    M_tot       = rhs.M_tot;                              // Accumulated number of phages
    dt          = rhs.dt;                                 // [hour]   Default size of the time step
    dT_c        = rhs.dT_c;                               // [hour]   Size of the time step which updates the phages
    dT          = rhs.dT;                                 // [hour]   Size of the time step which updates everything else

    dR          = rhs.dR;                                 // Depth of the infected layer (for spawning)
    nSamp       = rhs.nSamp;                              // Number of samples to save
    dGrid       = rhs.dGrid;                              // Spacing of grid (in units of critical radius R)
    dShell      = rhs.dShell;                             // [µm]     Thickness of agent layer (for biomass)

    n_0         = rhs.n_0;                                // [1/µm^3] Initial concentration of nutrient
    g_max       = rhs.g_max;                              // [1/hour] Maximal growth rate for the cells
    R           = rhs.R;                                  // [µm]     The length scale for division (Typical volume 1.3 µm^3)
    k           = rhs.k;                                  // [N*m]    Parameter for repulsive potential
    K           = rhs.K;                                  //          Michaels-Menten FACTOR for Monod growth
    gamma       = rhs.gamma;                              //          Probability to infect cell
    alpha       = rhs.alpha;                              //          Probability for phage to go lysogenic
    beta        = rhs.beta;                               //          Multiplication factor phage
    delta       = rhs.delta;                              // [1/hour] Rate of phage decay
    epsilon     = rhs.epsilon;                            //          Probability for offspring to turn resistant
    r           = rhs.r;                                  //          Constant used in the time-delay mechanism
    T_i         = rhs.T_i;                                // [hours]  Time when the phage infections begins (less than 0 disables phage infection)
    L_box       = rhs.L_box;                              // [µm]     Length of boundary condition box

    numB        = rhs.numB;                               // Current tally of suceptible cells
    numL        = rhs.numL;                               // Current tally of lysogenic cells
    numI        = rhs.numI;                               // Current tally of infected cells
    numD        = rhs.numD;                               // Current tally of dead cells
    numR        = rhs.numR;                               // Current tally of resistant cells

    eta         = rhs.eta;                                // Amount of division noise (width of gaussian)
    nu          = rhs.nu;                                 // Amount of displacement noise during division. x = x0 + rand(nu) etc.

    D_B         = rhs.D_B;                                // [µm^2/hour] Diffusion constant for the cells
    D_P         = rhs.D_P;                                // [µm^2/hour] Diffusion constant for the phage
    D_n         = rhs.D_n;                                // [µm^2/hour] Diffusion constant for the nutrient

    h_agar      = rhs.h_agar;                             // [µm]    Total height of the soft agar layers
    h_cell      = rhs.h_cell;                             // [µm]    Distance of the cell colony to the hard agar

    p           = rhs.p;                                  // Helper variable to control phage spawning

    debug       = rhs.debug;                              // The amount of information to print to terminal
    debugBool   = rhs.debugBool;
    exit        = rhs.exit;                               // Boolean to control early exit
    Time        = rhs.Time;                               // Counter for how many timesteps have passed
    maxStep     = rhs.maxStep;                            // Maximum time-steps to "step" when finegrained bool is false
    sigma       = rhs.sigma;                              // Maximum number of standard deviations to move in diffusions term

    nutrientField = rhs.nutrientField;                    // Boolean to toggle between nutrient field and just nutrients

    exportAny            = rhs.exportAny;                 //
    exportCellData       = rhs.exportCellData;            //
    exportColonySize     = rhs.exportColonySize;          // Booleans to control the export output
    exportPhageData      = rhs.exportPhageData;           //
    exportNutrient       = rhs.exportNutrient;            //

    singleInfectedCell   = rhs.singleInfectedCell;        // Boolean to enable phages via a single infected cell
    planarPhageInvasion  = rhs.planarPhageInvasion;       // Boolean to enable phages invading in a plane
    uniformPhageInvasion = rhs.uniformPhageInvasion;      // Boolean to enable phages invading uniformly from entire space
    manyInfectedCells    = rhs.manyInfectedCells;         // Boolean to enable surface infection and uniformly from entire spacefalse
    forcedThickness      = rhs.forcedThickness;           // "Boolean" to force an infected surface of a given number

    invasionType         = rhs.invasionType;              // Number to remember the invasion type

    reflectingBoundary   = rhs.reflectingBoundary;        // Boolean to enable the simple box boundary condition with reflective sides
    absorbingBoundary    = rhs.absorbingBoundary;         // Boolean to enable the simple box boundary condition with absorbing sides
    experimentBoundary   = rhs.experimentBoundary;        // Boolean to enable the advanced petri dish boundary conditions

    phageType  = rhs.phageType;                           // Contains the type of phage chosen (sets parameters according to De Paepe)
    rngSeed    = rhs.rngSeed;                             // The seed for the random number generator

    // Initilize the benchmarking vector
    for (int j = 0; j < 20; j ++) {
        benchMark[j] = rhs.benchMark[j];
    }

    res = rhs.res;

    // Copy random number generator
    rng = rhs.rng;

    // Copy cells, phages, nutrients and nutrient grid.
    cells           = rhs.cells;
    biomass         = rhs.biomass;
    phages          = rhs.phages;
    center          = rhs.center;
    nutrient        = rhs.nutrient;
    nutrient_grid   = rhs.nutrient_grid;
    CN              = rhs.CN;
    B               = rhs.B;
}


// Clean up /////////////////////////////////////////////////////////////////////////////
// Delete the data folder
void Simulation::DeleteFolder() {
    DeleteFolderTree(path.c_str());
}


// Delete folders recursively
void Simulation::DeleteFolderTree(const char* directory_name) {

    DIR*            dp;
    struct dirent*  ep;
    char            p_buf[512] = {0};


    dp = opendir(directory_name);
    if (dp == NULL) {
        return;
    }

    while ((ep = readdir(dp)) != NULL) {
        // Skip self dir "."
        if (strcmp(ep->d_name, ".") == 0 || strcmp(ep->d_name, "..") == 0) continue;

        sprintf(p_buf, "%s/%s", directory_name, ep->d_name);

        // Is the path a folder?
        struct stat s_buf;
        int IsDirectory = -1;
        if (stat(p_buf, &s_buf)){
            IsDirectory = 0;
        } else {
            IsDirectory = S_ISDIR(s_buf.st_mode);
        }

        // If it is a folder, go recursively into
        if (IsDirectory) {
            DeleteFolderTree(p_buf);
        } else {    // Else delete the file
            unlink(p_buf);
        }

    }

    closedir(dp);
    rmdir(directory_name);

}


// Destructor
Simulation::~Simulation() {

    // Close filestreams
    if (f_cells.is_open()) {
        f_cells.close();
    }
    if (f_colonySize.is_open()) {
        f_colonySize.close();
    }
    if (f_phages.is_open()) {
        f_phages.close();
    }
    if (f_n.is_open()) {
        f_n.close();
    }
    if (f_log.is_open()) {
        f_log.close();
    }
}