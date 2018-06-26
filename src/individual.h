#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

#include "global.h"
#include "fobjective.h"

class TIndividual
{
public:
	TIndividual();
	virtual ~TIndividual();

	void init_individual();
	void eval_objective();
	void norm_objective();

	bool operator<(const TIndividual &ind2);
	bool operator==(const TIndividual &ind2);
	void operator=(const TIndividual &ind2);

	void show_objective();
	void show_variable();

	vector<double> x_var;
	vector<double> y_obj;
	vector<double> n_obj; // normalized objectives
};

TIndividual::TIndividual()
{
	for (size_t i = 0; i < numVariables; i++)
		x_var.push_back(0.0);
	for (size_t i = 0; i < numObjectives; i++) {
		y_obj.push_back(0.0);
		n_obj.push_back(0.0);
	}
}

TIndividual::~TIndividual() {}

void TIndividual::init_individual()
{
	for (size_t i = 0; i < numVariables; i++)
		x_var[i] = lowBound + rnd_uni(&rnd_uni_init) * (uppBound - lowBound);
}

void TIndividual::eval_objective()
{
	objectives(x_var, y_obj);
}

void TIndividual::norm_objective()
{
	for (size_t i = 0; i < numObjectives; ++i) {
		if (fabs(nadirpoint[i] - idealpoint[i]) > 10e-10)
			n_obj[i] = (y_obj[i] - idealpoint[i]) / (nadirpoint[i] - idealpoint[i]);
		else
			n_obj[i] = (y_obj[i] - idealpoint[i]) / 10e-10;
	}
}

void TIndividual::operator=(const TIndividual &ind2)
{
    x_var = ind2.x_var;
	y_obj = ind2.y_obj;
	n_obj = ind2.n_obj;
}

bool TIndividual::operator<(const TIndividual &ind2)
{
	bool dominated = true;
	for (size_t i = 0; i < numObjectives; i++) {
		if (ind2.y_obj[i] < y_obj[i])
			return false;
	}
	if (ind2.y_obj == y_obj)
		return false;
	return dominated;
}

bool TIndividual::operator==(const TIndividual &ind2)
{
	if (ind2.y_obj == y_obj)
		return true;
	else
		return false;
}

void TIndividual::show_objective()
{
	for (size_t i = 0; i < numObjectives; i++)
		printf("%f ", y_obj[i]);
	printf("\n");
}

void TIndividual::show_variable() {
	for (size_t i = 0; i < numVariables; i++)
		printf("%f ", x_var[i]);
	printf("\n");
}

/* ******************************************************************
 * TParetofronts class
 * ******************************************************************/

class TParetofront
{
public:
	TParetofront() {};
	virtual ~TParetofront() {};

	typedef vector<int> TPFront;
	typedef vector<TPFront> TPfronts;
};

#endif /* INDIVIDUAL_H_ */
