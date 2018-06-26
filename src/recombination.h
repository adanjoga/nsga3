#ifndef _RECOMBINATION_H_
#define _RECOMBINATION_H_

#include "global.h"

/* Routine for real polynomial mutation of an T */
template<class T>
void realmutation(T &ind, double rate) {
	double rnd, delta1, delta2, mut_pow, deltaq;
	double y, yl, yu, val, xy;
	double eta_m = id_mu;

	for (int j = 0; j < numVariables; j++) {
		if (rnd_uni(&rnd_uni_init) <= rate) {
			y = ind.x_var[j];
			yl = lowBound;
			yu = uppBound;
			delta1 = (y - yl) / (yu - yl);
			delta2 = (yu - y) / (yu - yl);
			rnd = rnd_uni(&rnd_uni_init);
			mut_pow = 1.0 / (eta_m + 1.0);
			if (rnd <= 0.5) {
				xy = 1.0 - delta1;
				val = 2.0 * rnd + (1.0 - 2.0 * rnd) * (pow(xy, (eta_m + 1.0)));
				deltaq = pow(val, mut_pow) - 1.0;
			} else {
				xy = 1.0 - delta2;
				val = 2.0 * (1.0 - rnd)
						+ 2.0 * (rnd - 0.5) * (pow(xy, (eta_m + 1.0)));
				deltaq = 1.0 - (pow(val, mut_pow));
			}
			y = y + deltaq * (yu - yl);
			if (y < yl)
				y = yl;
			if (y > yu)
				y = yu;
			ind.x_var[j] = y;
		}
	}
	return;
}

/* Routine for real variable SBX crossover */
template<class T>
void realbinarycrossover(T &parent1, T &parent2, T& child1, T& child2) {
	double rand;
	double y1, y2, yl, yu;
	double c1, c2;
	double alpha, beta, betaq;
	double eta_c = id_cx;
	if (rnd_uni(&rnd_uni_init) <= 1.0) {
		for (int i = 0; i < numVariables; i++) {
			if (rnd_uni(&rnd_uni_init) <= 0.5) {
				if (fabs(parent1.x_var[i] - parent2.x_var[i]) > EPS) {
					if (parent1.x_var[i] < parent2.x_var[i]) {
						y1 = parent1.x_var[i];
						y2 = parent2.x_var[i];
					} else {
						y1 = parent2.x_var[i];
						y2 = parent1.x_var[i];
					}
					yl = lowBound;
					yu = uppBound;
					rand = rnd_uni(&rnd_uni_init);
					beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(eta_c + 1.0));
					if (rand <= (1.0 / alpha)) {
						betaq = pow((rand * alpha), (1.0 / (eta_c + 1.0)));
					} else {
						betaq = pow((1.0 / (2.0 - rand * alpha)),
								(1.0 / (eta_c + 1.0)));
					}
					c1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
					beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
					alpha = 2.0 - pow(beta, -(eta_c + 1.0));
					if (rand <= (1.0 / alpha)) {
						betaq = pow((rand * alpha), (1.0 / (eta_c + 1.0)));
					} else {
						betaq = pow((1.0 / (2.0 - rand * alpha)),
								(1.0 / (eta_c + 1.0)));
					}
					c2 = 0.5 * ((y1 + y2) + betaq * (y2 - y1));
					if (c1 < yl)
						c1 = yl;
					if (c2 < yl)
						c2 = yl;
					if (c1 > yu)
						c1 = yu;
					if (c2 > yu)
						c2 = yu;
					if (rnd_uni(&rnd_uni_init) <= 0.5) {
						child1.x_var[i] = c2;
						child2.x_var[i] = c1;
					} else {
						child1.x_var[i] = c1;
						child2.x_var[i] = c2;
					}
				} else {
					child1.x_var[i] = parent1.x_var[i];
					child2.x_var[i] = parent2.x_var[i];
				}
			} else {
				child1.x_var[i] = parent1.x_var[i];
				child2.x_var[i] = parent2.x_var[i];
			}
		}
	} else {
		for (int i = 0; i < numVariables; i++) {
			child1.x_var[i] = parent1.x_var[i];
			child2.x_var[i] = parent2.x_var[i];
		}
	}
	return;
}

#endif
