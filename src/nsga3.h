#ifndef NSGA3_H_
#define NSGA3_H_

#include "global.h"
#include "individual.h"
#include "refpoint.h"
#include "recombination.h"
#include "common.h"

class TNSGA3 {

public:
	TNSGA3();
	virtual ~TNSGA3();

	int popsize;					// population size
	vector<TIndividual> parent_pop;	// parent population
	vector<TIndividual> child_pop;	// child population
	vector<TIndividual> mixed_pop;	// combined population
	vector<TRefpoint> refepoints;	// reference points

	void init_population();			 	// initialize the population
	void init_refpoints2D(); 			// initialize the reference points
	void init_refpointsND();			// generate the reference points
	void gen_refpointsND(vector<TRefpoint> *rps, TRefpoint *rp, size_t left, size_t count);
	void merge_populations();			// combination of parent and child populations

	void update_idealpoint(TIndividual &indv);
	void run(size_t ngen, size_t run);
	void evolution();
	TParetofront::TPfronts nondominatedsort();

	void normalize(TParetofront::TPfronts Fronts);
	void associate(TParetofront::TPfronts Fronts);
	void niching();

	void save_front(char file_name[512]);
};

TNSGA3::TNSGA3() {
	popsize = 0;
	idealpoint = new double[numObjectives];
	nadirpoint = new double[numObjectives];
	for (int i = 0; i < numObjectives; i++) {
		idealpoint[i] = MAX_DOUBLE;
		nadirpoint[i] = MIN_DOUBLE;
	}
}

TNSGA3::~TNSGA3() {
	delete[] idealpoint;
}

// initialize a set of evenly-distributed weight vectors for 2 objectives
void TNSGA3::init_refpoints2D()
{
	for (size_t i = 0; i <= nrefpoints_2D; i++) {
		TRefpoint rpoint;
		vector<size_t> array;
		array.push_back(i);
		array.push_back(nrefpoints_2D - i);
		for (size_t j = 0; j < array.size(); j++)
			rpoint.refpoint.push_back(1.0 * array[j] / nrefpoints_2D);
		refepoints.push_back(rpoint);
	}
}
// initialize a set of evenly-distributed weight vectors for <2 objectives
 void TNSGA3::init_refpointsND()
 {
	TRefpoint rpoint(numObjectives);
	gen_refpointsND(&refepoints, &rpoint, p_boundary, 0);

	if (p_inside > 0) {
		vector<TRefpoint> inside_rps;	// inside reference points
		gen_refpointsND(&refepoints, &rpoint, p_inside, 0);

		double center = 1.0/numObjectives;
		for (size_t i = 0; i < inside_rps.size(); ++i) {
			for (size_t j = 0; j < numObjectives; ++j)
				inside_rps[i].refpoint[j] = (center + inside_rps[i].refpoint[j])/2;
			refepoints.push_back(inside_rps[i]);
		}
	}
 }

// generate a set of evenly-distributed weight vectors for <2 objectives
void TNSGA3::gen_refpointsND(vector<TRefpoint> *rps, TRefpoint *rp, size_t left, size_t count)
{
	if (count == numObjectives-1) {
		rp->refpoint[count] = 1.0 * left / p_boundary;
		rps->push_back(*rp);
	}
	else {
		for (size_t i = 0; i <= left; ++i) {
			rp->refpoint[count] = 1.0 * i / p_boundary;
			gen_refpointsND(rps, rp, left-i, count+1);
		}
	}
}

void TNSGA3::init_population()
{
	for (size_t i = 0; i < popsize; ++i) {
		TIndividual indv;
		indv.init_individual();
		indv.eval_objective();
		parent_pop.push_back(indv);
	}
}

void TNSGA3::merge_populations()
{
	for (size_t i = 0; i < popsize; ++i)
		mixed_pop.push_back(parent_pop[i]);

	for (size_t i = 0, j = popsize; i < popsize; ++i, j++)
		mixed_pop.push_back(child_pop[j]);

}

TParetofront::TPfronts TNSGA3::nondominatedsort()
{
	int rank = 1, count = 0;
	vector<int> ranks(mixed_pop.size(), 0);
	TParetofront::TPfronts fronts;

	while (count < mixed_pop.size()) {
		TParetofront::TPFront cfront;
		for (size_t i = 0; i < mixed_pop.size(); ++i) {
			if (ranks[i] > 0)
				continue;
			bool isdominated = false;
			for (size_t j = 0; j < cfront.size(); ++j) {
				if (mixed_pop[cfront[j]] < mixed_pop[i]) {
					isdominated = true;
					break;
				} else if (mixed_pop[i] < mixed_pop[cfront[j]]) {
					cfront.erase(cfront.begin() + j);
					j = j - 1;
				}
			}
			if (!isdominated) cfront.push_back(i);
		}
		for (size_t i = 0; i < cfront.size(); ++i) ranks[cfront[i]] = rank;
		fronts.push_back(cfront);
		count += cfront.size();
		rank++;
	}

	// identify the last front index (fidx)
	int fidx = 0, block = 0;
	while (block < popsize) block += fronts[fidx++].size();

	// remove useless fronts
	fronts.erase(fronts.begin() + fidx, fronts.end());

	return fronts;
}

void TNSGA3::normalize(TParetofront::TPfronts Fronts) {
	// compute the ideal point
	for (size_t i = 0; i < Fronts[0].size(); ++i) {
		size_t idx = Fronts[0][i];
		for (size_t j = 0; j < numObjectives; ++j) {
			if (mixed_pop[idx].y_obj[j] < idealpoint[j])
				idealpoint[j] = mixed_pop[idx].y_obj[j];
		}
	}
	// compute the Nadir point
	for (size_t i = 0; i < Fronts.size(); ++i) {
		for (size_t j = 0; j < Fronts[i].size(); ++j) {
			size_t idx = Fronts[i][j];
			for (size_t k = 0; k < numObjectives; ++k) {
				if (mixed_pop[idx].y_obj[k] > nadirpoint[k])
					nadirpoint[k] = mixed_pop[idx].y_obj[k];
			}
		}
	}
	// Normalize objectives
	for (size_t i = 0; i < Fronts.size(); ++i)
		for (size_t j = 0; j < Fronts[i].size(); ++j)
			mixed_pop[Fronts[i][j]].norm_objective();
}

void TNSGA3::associate(TParetofront::TPfronts Fronts)
{
	for (size_t i = 0; i < refepoints.size(); ++i)
		refepoints[i].clear_refpoint();

	for (size_t i = 0; i < Fronts.size(); ++i) {
		for (size_t j = 0; j < Fronts[i].size(); ++j) {

			double mindist = MAX_DOUBLE;
			size_t idx = refepoints.size();
			for (size_t r = 0; r < refepoints.size(); ++r) {
				double pdist = perpendicular_distance(refepoints[r].refpoint,
						mixed_pop[Fronts[i][j]].n_obj);
				if (pdist < mindist) {
					mindist = pdist;
					idx = r;
				}
			}
			if (i + 1 == Fronts.size())
				refepoints[idx].add_member(Fronts[i][j], mindist);
			else
				refepoints[idx].sum_member();
		}
	}
}

void TNSGA3::niching()
{
	vector<bool> rp_isable (refepoints.size(), true);
	while (parent_pop.size() < popsize) {

		size_t rp_idx = find_refpoint(refepoints, rp_isable);
		int ind_idx = select_member(refepoints[rp_idx]);

		if (ind_idx < 0) rp_isable[rp_idx] = false;
		else {
			refepoints[rp_idx].sum_member();
			refepoints[rp_idx].remove_member(ind_idx);
			parent_pop.push_back(mixed_pop[ind_idx]);
		}
	}
}

void TNSGA3::evolution()
{
	child_pop.clear(); // clear the child population
	for (size_t i = 0; i < popsize; ++i) {
		int p1 = int(popsize * rnd_uni(&rnd_uni_init));
		int p2 = int(popsize * rnd_uni(&rnd_uni_init));

		TIndividual child1, child2;
		realbinarycrossover(parent_pop[p1], parent_pop[p2], child1, child2);
		realmutation(child1, 1.0 / numVariables);
		realmutation(child2, 1.0 / numVariables);

		child1.eval_objective();
		child2.eval_objective();

		child_pop.push_back(child1);
		child_pop.push_back(child2);
	}

	// mixed_pop: combination of parents and child population
	mixed_pop.clear();
	merge_populations();

	// non-dominated sorting algorithm
	TParetofront::TPfronts Fronts = nondominatedsort();

	// copy the first F-1 fronts to the next population
	parent_pop.clear();
	for (int i = 0; i < Fronts.size() - 1; ++i)
		for (int j = 0; j < Fronts[i].size(); ++j)
			parent_pop.push_back(mixed_pop[Fronts[i][j]]);

	// finish if the population is complete
	if (parent_pop.size() == popsize) return;

	// normalization of objectives (algorithm 2:normalize function)
	normalize(Fronts);

	// association of members (algorithm 3: associate function)
	associate(Fronts);

	// selection of k individual (algorithm 4: niching function)
	niching();
}

void TNSGA3::run(size_t ngen, size_t run) {

	if (numObjectives == 2) init_refpoints2D();
	else init_refpointsND();

	popsize = refepoints.size();
	while (popsize % 4) popsize++;	// population size must be multiple of 4

	init_population();
	for (size_t g = 0; g < ngen; ++g) evolution();

	char file_name[512];
	sprintf(file_name,"ParetoFront/NSGA3_%s_%dD_R%lu.data", strTestInstance, numObjectives, run);
	save_front(file_name);

	parent_pop.clear();
}

void TNSGA3::save_front(char file_name[512])
{
	std::fstream fout;
	fout.open(file_name, std::ios::out);
	for (size_t i = 0; i < parent_pop.size(); i++) {
		for (int j = 0; j < numObjectives; j++)
			fout << parent_pop[i].y_obj[j] << "  ";
		fout << "\n";
	}
	fout.close();
}


#endif /* NSGA3_H_ */
