#include "Simulation.hpp"

using namespace std;

int main(int argc, char** argv){


	double gamma[9] = {  pow(10,6), pow(10,5.5), pow(10,5), pow(10,4.5), pow(10,4), pow(10,3.5), pow(10,3), pow(10,2.5), pow(10,2)};

	double L = 25;


	#pragma omp parallel for
	for (int i = 0; i < 9; i++) {


		// Load simulation module
		Simulation s(1);

		s.Quiet();

		s.SetFolderNumber(i);

		// Set the invasion time
		s.PhageInvasionStartTime(0);

		// Set the invasion type
		s.PhageInvasionType(3);			// 1: Single Infected cell, 2: Planar Phage Invasion, 3: Uniform Phage Invasion, 4: Many Infected Cells
        s.BoundaryType(2);

        // Sets initial density of the phage
        s.PhageInitialDensity(1e4 / pow(L, 3));

        // Sets initial density of the bacteria
		s.CellInitialDensity(1/pow(L,2));

		// Set the infection time
		s.PhageInfectionRate(0.0);
        s.PhageDecayRate(0.0);
		s.MaxGrowthRate(0.0);

        s.SetSamples(1000);

		// s.AdsorptionParameter(1);
		// s.AdsorptionParameter(0.01);
		s.PhageAdsorptionParameter(gamma[i]);

		// Set the data to export
		// s.ExportCellData();
		s.ExportColonySize();

		// Start the simulation
		s.Run(10);

		}

	return 0;
}