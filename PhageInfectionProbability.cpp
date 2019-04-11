#include "Simulation.hpp"

using namespace std;

int main(int argc, char** argv){

    double gamma = 1e6;

    for (int i = 0; i < argc; i++) {
        if (i == 1) {
            gamma = atof(argv[i]);
        }
    }

    // Set the sizes the colonies should reach
    vector<int> ColonySize = {10, 32, 100, 316, 1000};

    // Set the percentage of infected bacteria
    vector<double> f = {0.2, 0.4, 0.6, 0.8};

    double L = 25;

    // #pragma omp parallel for
    for (int i = 0; i < ColonySize.size(); i++) {

        // Load master simulation
        Simulation m(1);
        m.Quiet();

        // Set a random seed
        m.SetRngSeed(0);

        // Sets initial density of the bacteria
        m.CellInitialDensity(1 / pow(L, 2));

        // Set phage properties
        string phage = "P1vir";
        m.PhageType(phage);

        m.PhageAdsorptionParameter(gamma);

        // Set boundary conditions
        m.BoundaryType(2); // Reflecting

        // Set maximum colony size:
        m.MaxColonySize(ColonySize[i]);

        // Run until size is reached
        m.Run(1e-5);
        double T = 1e-5;
        while (m.NumberOfUninfectedCells() < ColonySize[i]) {
            int err = m.Run(1e-5);
            T += 1e-5;

            if (err > 0) {
                break;
            }
            if (T > 10) {
                break;
            }
        }

        cout << "ColonySize was set to: " << ColonySize[i] << ".\t size obtained is: " << m.NumberOfUninfectedCells() << endl;
        cout << "gamma was set to: " << gamma << ".\t gamma obtained is: " << m.GetGamma() << endl;

        // #pragma omp parallel for
        for (int j = 0; j < f.size(); j++) {

            // Copy the simulation state
            Simulation c(m);
            c.ClearErrors();

            c.TimeStepSkip(1);

            // Set the number of infected in surface layer
            c.PhageInvasionType(6, (int)round(f[j] * ColonySize[i]));

            // Set the infection time to be next timestep
            double T_i = static_cast<double>(c.GetTime() + 1) * 1e-6;
            c.PhageInvasionStartTime(T_i);

            // Stop colony growth
            c.MaxGrowthRate(0.0);

            // Run time forward by dT to check spawning
            c.Run(1e-5);
            cout << "forcedInfected was set to: " << (int)round(f[j] * ColonySize[i]) << ".\t size obtained is: " << m.NumberOfUninfectedCells() - c.NumberOfUninfectedCells() << endl;

            // Spawnphages to diffuse around colony
            int N_hit = 0;
            int N_try = 0;
            int k_max = 5e4;

            for (int k = 0; k < k_max; k++) {

                N_try++;

                // Copy the simulation state
                Simulation t(c);

                // Set uniform phage invasion
                t.PhageInvasionType(3);

                // Set the infection time to be next timestep
                double T_i = static_cast<double>(t.GetTime()+1)*1e-6;
                t.PhageInvasionStartTime(T_i);

                // Set density so we get one phage
                t.PhageInitialDensity(1/pow(L,3));

                t.PhageAdsorptionParameter(gamma);

                // Set the rngseed
                t.SetRngSeed(k);

                // Run time forward by dT to check spawning
                t.Run(1e-5);

                int N_it = 0;
                while (t.NumberOfPhages() > 0) {
                    t.Run(1e-2);
                    N_it++;
                }

                if (t.NumberOfUninfectedCells() < c.NumberOfUninfectedCells()) {
                    N_hit++;
                }

                if (N_hit >= 1000) {
                    k = k_max;
                }
            }

            cout << "\t# of new infections: " << N_hit << " out of "<< N_try << " attempts" << endl;
            cout << "\tThe probability of hitting a non-infected cell is: " << (double)N_hit/(double)N_try << endl << endl;

        }
        cout << endl;
        cout << endl;
	}

	return 0;
}
