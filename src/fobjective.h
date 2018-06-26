#ifndef FOBJECTIVE_H_
#define FOBJECTIVE_H_

#include "global.h"

#define pi   3.1415926
#define SQR2  sqrt(2)

void objectives(vector<double> x_var, vector<double> &y_obj)
{

	if(!strcmp(strTestInstance,"DTLZ1_temp"))
	{
		double g = 0;
		for (int n = 2; n < numVariables; n++)
			g = g + pow(x_var[n] - 0.5, 2) - cos(20 * pi * (x_var[n] - 0.5));
		g = 100 * (numVariables - 2 + g);
		y_obj[0] = (1 + g) * x_var[0] * x_var[1];
		y_obj[1] = (1 + g) * x_var[0] * (1 - x_var[1]);
		y_obj[2] = (1 + g) * (1 - x_var[0]);
	}

	if(!strcmp(strTestInstance,"DTLZ1"))
	{
		double k = numVariables - numObjectives + 1;
		double g = 0.0;

		for (int i = numVariables - k; i < numVariables; i++)
			g += (x_var[i] - 0.5) * (x_var[i] - 0.5)
					- cos(20.0 * pi * (x_var[i] - 0.5));

		g = 100.0 * (k + g);
		for (int i = 0; i < numObjectives; i++)
			y_obj[i] = (1.0 + g) * 0.5;

		for (int i = 0; i < numObjectives; i++) {
			for (int j = 0; j < numObjectives - (i + 1); j++)
				y_obj[i] *= x_var[j];
			if (i != 0) {
				int aux = numObjectives - (i + 1);
				y_obj[i] *= 1 - x_var[aux];
			}
		}
	}

	if (!strcmp(strTestInstance, "DTLZ2"))
	{
		double k = numVariables - numObjectives + 1;
		double g = 0.0;

		for (int i = numVariables - k; i < numVariables; i++)
			g += (x_var[i] - 0.5) * (x_var[i] - 0.5);

		for (int i = 0; i < numObjectives; i++)
			y_obj[i] = 1.0 + g;

		for (int i = 0; i < numObjectives; i++) {
			for (int j = 0; j < numObjectives - (i + 1); j++)
				y_obj[i] *= cos(x_var[j] * 0.5 * pi);
			if (i != 0) {
				int aux = numObjectives - (i + 1);
				y_obj[i] *= sin(x_var[aux] * 0.5 * pi);
			}
		}
	}

	if (!strcmp(strTestInstance, "DTLZ3"))
	{
		double k = numVariables - numObjectives + 1;
		double g = 0.0;

		for (int i = numVariables - k; i < numVariables; i++)
			g += (x_var[i] - 0.5) * (x_var[i] - 0.5)
					- cos(20.0 * pi * (x_var[i] - 0.5));

		g = 100.0 * (k + g);
		for (int i = 0; i < numObjectives; i++)
			y_obj[i] = 1.0 + g;

		for (int i = 0; i < numObjectives; i++) {
			for (int j = 0; j < numObjectives - (i + 1); j++)
				y_obj[i] *= cos(x_var[j] * 0.5 * pi);
			if (i != 0) {
				int aux = numObjectives - (i + 1);
				y_obj[i] *= sin(x_var[aux] * 0.5 * pi);
			}
		}
	}

	if (!strcmp(strTestInstance, "DTLZ4"))
	{
		double k = numVariables - numObjectives + 1;
		double g = 0.0, alpha = 100.0;

		for (int i = numVariables - k; i < numVariables; i++)
			g += (x_var[i] - 0.5) * (x_var[i] - 0.5);

		for (int i = 0; i < numObjectives; i++)
			y_obj[i] = 1.0 + g;

		for (int i = 0; i < numObjectives; i++) {
			for (int j = 0; j < numObjectives - (i + 1); j++)
				y_obj[i] *= cos(pow(x_var[j], alpha) * (pi / 2.0));
			if (i != 0) {
				int aux = numObjectives - (i + 1);
				y_obj[i] *= sin(pow(x_var[aux], alpha) * (pi / 2.0));
			}
		}
	}
}

#endif /* FOBJECTIVE_H_ */
