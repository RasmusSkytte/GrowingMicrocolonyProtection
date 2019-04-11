#include "Simulation.hpp"

using namespace std;

int main(int argc, char** argv){


	double D_P[5] = {1000,  2000,  4000,  9000, 17000 };
	double T_i[5] = {2.00,  2.00,  2.00,  2.00, 2.00  };

	double L = 65;
	double P0 = 0.01;

	// Searching steps
	double dT[3] = { 0.5, 0.1, 0.01 };

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

		m.TimeStepSkip(33);

		// Sets initial density of the bacteria
		m.CellInitialDensity(1/pow(L,2));

		// Set phage properties
		string phage = "P1vir";
		m.PhageType(phage);
		m.PhageDiffusionConstant(D_P[i]);

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
				stream << "Initial time too small! (T_i[" << i << "] = " << T_i[i] << ")" << endl;
				cout << stream.str();
				break;
			}

			// Set time-step size
			double dt = dT[it];

			// Break for loop if convergence is found
			if ((err == 0) && (it == 3)) {
				break;
			}

			// Safety break, if things go wrong..
			if ((k > 10) or (l > 10)) {
				break;
			}

			// Make copy for this resolution
			Simulation s(c);
			s.SetFolderNumber(1000*i + 100*j + 10*k + l + 20000);

			// Set phage properties
			s.PhageInvasionStartTime(T_i[i]);
			s.PhageInvasionType(3);			// 1: Single Infected cell, 2: Planar Phage Invasion, 3: Uniform Phage Invasion, 4: Many Infected Cells
			s.PhageInitialDensity(P0);

			s.BoundaryType(2); 				// 1: Absorbing, 2: Reflecting, 3: Experimental

			// Set the data to export
			s.ExportCellData();
			s.ExportColonySize();
			s.ExportNutrient();

			// Run for 1 hours to get inital number of uninfected cells
			// 1 hour, since the first wave of infected cells should lyse after 40 minutes.
			err = s.Run(1);

			// Get initial cell count
			int N0 = s.NumberOfUninfectedCells();

			// Run simulation for 2 additional hours
			for (int j = 0; j < 4; j++) {

				// Exit if flag is 1
				if (err == 1) {
					break;
				}

				// Exit if N0 = 0
				if (N0 == 0) {
					break;
				}

				// Run simulation for 0.1 hours
				err = s.Run(0.5);

				// Exit if there are fewer cells than initally:
				if (N0 > s.NumberOfUninfectedCells()) {
					err = 1;
					break;
				}

				// Exit if growth is sufficiently fast:
				if (1.65*N0 < s.NumberOfUninfectedCells()) {
					err = 0;
					break;
				}

				// Update N0
				N0 = s.NumberOfUninfectedCells();	// Get new value

			}

			// err is now either 1 or 0
			if (err == 1) { // The run was unsuccesful, the time should be driven forward

				m = c;		// Merge the copy state with the master state
				c.Run(dt);	// Run the copy forward in time
				T_i[i] += dt;

				// Increment counter
				if (it == 0) {
					j++;
				} else if (it == 1) {
					k++;
				} else if (it == 2) {
					l++;
				}

				// Remove the data folder
				s.DeleteFolder();

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
					m.Run(dT[it]);
					T_i[i] -= dt;		// Subtract the test increment
					T_i[i] += dT[it];	// Add the next test increment

					// Load the master into the copy
					c = m;

					// Remove the data folder
					s.DeleteFolder();

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
