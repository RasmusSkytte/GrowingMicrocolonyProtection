#include "Simulation.hpp"

using namespace std;

int main(int argc, char** argv){

	// Infection rates
	double r[5]   = { 0.4,  0.5,  0.6,  0.7,  0.8 };
	double T_i[5] = { 0.0,  0.0,  0.0,  0.0,  0.0 };

	double L = 65;
	double b = 0;
	double P0 = 0.001;

	// Searching steps
	double dT[3] = { 0.25, 0.05, 0.01 };

	int j = 0;

	if (j == 0) { 	// No Lysogeny (b = 100, P0 = 0.01)
		T_i[0]  = 1.0;
		T_i[1]  = 1.0;
		T_i[2]  = 1.0;
		T_i[3]  = 1.0;
		T_i[4]  = 1.0;

		b = 100;
	}
	if (j == 1) { 	// No Lysogeny (b = 1000, P0 = 0.01)
		T_i[0]  = 1.0;
		T_i[1]  = 1.0;
		T_i[2]  = 1.0;
		T_i[3]  = 1.0;
		T_i[4]  = 1.0;

		b = 1000;
	}
	if (j == 2) { 	// No Lysogeny (b = 400, P0 = 0.01)
		T_i[0]  = 1.0;
		T_i[1]  = 1.0;
		T_i[2]  = 1.0;
		T_i[3]  = 1.0;
		T_i[4]  = 1.0;

		b = 400;
	}

	#pragma omp parallel for
	for (int i = 0; i < 5; i++) {

		// Start timer
		time_t timer;
		time(&timer);

		// Load master simulation
		Simulation m(1);
		m.Quiet();

		// Set a random seed
		m.SetRngSeed(0);

		m.TimeStepSkip(50);
		m.MaxStep(2);	// Maximum step of 2 sigma

		// Sets initial density of the bacteria
		m.CellInitialDensity(1/pow(L,2));

		// Set phage properties
		string phage = "P1vir";
		m.PhageType(phage);
		m.PhageInfectionRate(r[i]);
		m.PhageBurstSize(b);

		// Run until start time.
		m.Run(T_i[i]);

		// Copy the simulation state
		Simulation c(m);

		// Start iteration algorithm
		int it = 0;
		int j = 0;
		int k = 0;
		int l = 0;
		int err = 1;
		while (j < 10) {

			// Verify j,k,l is larger than 1
			if ((j < 0) or (k < 0) or (l < 0)) {
				stringstream stream;
				stream << "Initial time large small! (T_i[" << i << "] = " << T_i[i] << ")" << endl;
				cout << stream.str();
				break;
			}

			// Safety break, if things go wrong..
			if ((k > 10) or (l > 10)) {
				break;
			}

			// Set time-step size
			double dt = dT[it];

			// Break for loop if convergence is found
			if ((err == 0) && (it == 3)) {
				break;
			}

			// Make copy for this resolution
			Simulation s(c);
			s.SetFolderNumber(1000*i + 100*j + 10*k + l + 10000);

			// Set phage properties
			s.PhageInvasionStartTime(T_i[i]);
			s.PhageInvasionType(5);			// 1: Single Infected cell, 2: Planar Phage Invasion, 3: Uniform Phage Invasion, 4: Many Infected Cells, 5: Single Infected cell + Uniform Phage Invasion
			s.PhageInitialDensity(P0);

			s.BoundaryType(1); 				// 1: Absorbing, 2: Reflecting, 3: Experimental

			// Set the data to export
			s.ExportCellData();
			s.ExportColonySize();
			s.ExportNutrient();

			// Run for 1 hours to get inital number of uninfected cells
			// 1 hour, since the first wave of infected cells should lyse after 40 minutes.
			err = s.Run(1);

			// // Get initial cell count
			// int N0 = s.NumberOfUninfectedCells();

			// Get initial volume
			double V = s.Biomass();

			// Run simulation for 4 additional hours
			for (int n = 0; n < 8; n++) {

				// Exit if flag is 1
				if (err == 1) {
					break;
				}

				// // Exit if N0 = 0
				// if (N0 == 0) {
				// 	break;
				// }
				// Exit if there are no more cells
				if (V < 0.1) {
					break;
				}

				// Run simulation for 0.1 hours
				err = s.Run(0.5);

				// Exit if flag is 1
				if (err == 1) {
					break;
				}

				// // Exit if there are fewer cells than initally:
				// if (N0 > s.NumberOfUninfectedCells()) {
				// 	err = 1;
				// 	break;
				// }
				if (V > s.Biomass()) {
					err = 1;
					break;
				}

				// // Exit if growth is sufficiently fast:
				// if (1.65*N0 < s.NumberOfUninfectedCells()) {
				// 	err = 0;
				// 	break;
				// }

				// // Update N0
				// N0 = s.NumberOfUninfectedCells();	// Get new value

				V = s.Biomass();

			}

			// err is now either 1 or 0
			if (err == 1) { // The run was unsuccesful, the time should be driven forward

				m = c;		// Merge the copy state with the master state
				T_i[i] += dt;
				c.PhageInvasionStartTime(T_i[i]);
				c.Run(0);

				// Increment counter
				if (it == 0) {
					j++;
				} else if (it == 1) {
					k++;
				} else if (it == 2) {
					l++;
				}

				// Remove the data folder
				// s.DeleteFolder();

			} else { // The run was succesful, and the time step can be decreased and copy can be loaded

				// Decrement counter
				if (it == 0) {
					j--;
					k++;
				} else if (it == 1) {
					k--;
					l++;
				} else if (it == 2) {
					l--;
				}

				// Increase iterator
				it++;

				// If more iterations are to follow
				if (it < 3) {

					// Reset err
					err = 1;

					// Run the master simulation forward in time by dT
					T_i[i] -= dt;		// Subtract the test increment
					T_i[i] += dT[it];	// Add the next test increment
					m.PhageInvasionStartTime(T_i[i]);
					m.Run(0);

					// Load the master into the copy
					c = m;

					// Remove the data folder
					// s.DeleteFolder();

				} else {

					// Stop timer, and write result to datafolder
					float seconds = difftime(time(NULL),timer);
					float hours   = floor(seconds/3600);
					float minutes = floor(seconds/60);
					minutes -= hours*60;
					seconds -= minutes*60 + hours*3600;

					std::ofstream f_out;
					f_out.open(s.GetPath() + "/Completed.txt",fstream::out);
					f_out << "\tSimulation complete after ";
					if (hours > 0.0)   f_out << hours   << " hours and ";
					if (minutes > 0.0) f_out << minutes << " minutes and ";
					f_out  << seconds << " seconds." << endl;
					f_out.close();

					f_out.open(s.GetPath() + " - Completed.txt",fstream::out);
					f_out.close();
				}
			}
		}
	}


	return 0;
}
