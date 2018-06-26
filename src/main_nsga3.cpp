/*============================================================================
 // Name        : nsga3.cpp
 // Author      : Adán José-García
 // Version     : v1.0
 // Description : http://dx.doi.org/10.1109/TEVC.2013.2281535
 // Test the nsga3 algorithm and DTLZ2 problem varying the number of objectives
 //============================================================================*/

#include "global.h"
#include "nsga3.h"

int main() {

	size_t num_exp 	= 20;	// total number of experiments (runs)
	numVariables 	= 12; 	// number of variables (M+k-1) = (3+10-1), k=10
	strcpy(strTestInstance,"DTLZ2");

	ifstream indata("TestDTLZ2.txt");
	if (!indata) {
		cerr << "Error: file could not be opened" << endl;
		exit(1);
	}

	char temp[1024];
	while(!indata.eof()) {
		indata >> temp >> numObjectives;
		indata >> temp >> max_gen;
		indata >> temp >> p_boundary;
		indata >> temp >> p_inside;

		for (size_t run = 1; run <= num_exp; ++run) {
			printf(" Running experiment %lu for %s problem with %d objectives\n",
					run, strTestInstance, numObjectives);

			seed = (seed + 111) % 1235;
			rnd_uni_init = -(long) seed;

			TNSGA3 NSGA3;
			NSGA3.run(max_gen, run);
		}
	}
	indata.close();
	return 1;
}
