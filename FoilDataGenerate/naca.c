#include <stdio.h>
#include <math.h>
#include "draw_foil.hpp"
#define m 0.2025
#define p 0.15
#define k1 15.957
#define t 0.12
#define N 30
#define PI 3.1415926
typedef struct Point Point;
struct Point {
	double x;
	double y;
};

double y_t[N + 1] = { 0.0 };
double y_c[N + 1] = { 0.0 };
double x_u[N + 1], y_u[N + 1], x_l[N + 1], y_l[N + 1];
//Point point_upper[N], point_lower[N];

void generate_foil() {
	double x[N + 1] = { 0.0 };
	double theta[N + 1] = { 0.0 };	//theta is defined as arctan(partial y_c / partial x)
	for (int i = 0; i < N; i++) {
		x[i] = 1.0 * i / N;		//nondimensionalize x
		x[i] = pow(x[i], 2.5) * exp(x[i] * x[i] - 1);
		printf("x[%d] = %lf\n", i, x[i]);
		y_t[i] = t / 0.2 * (0.2969 * sqrt(x[i]) - 0.126 * x[i] - 0.3516 * x[i] * x[i] +
			0.2843 * pow(x[i], 3) - 0.1015 * pow(x[i], 4));
		if (x[i] <= p) {
			y_c[i] = k1 / 6 * (pow(x[i], 3) - 3 * m * pow(x[i], 2) + m * m * (3 - m) * x[i]);
			theta[i] = atan( k1 / 6 * (3 * pow(x[i], 2) - 6 * m * x[i] + m * m * (3 - m)));
		}
		else {
			y_c[i] = k1 / 6 * pow(m, 3) * (1 - x[i]);
			theta[i] = atan(-k1 / 6 * pow(m, 3));
		}
		x_u[i] = x[i] - y_t[i] * sin(theta[i]);
		y_u[i] = y_c[i] + y_t[i] * cos(theta[i]);
		x_l[i] = x[i] +  y_t[i] * sin(theta[i]);
		y_l[i] = y_c[i] - y_t[i] * cos(theta[i]);
	}
	x_u[N] = 1.0; x_l[N] = 1.0;
	y_u[N] = 0.0; y_l[N] = 0.0;
}
/* s = 'u' means to output upper foil data
while 'l' means lower;
format in each line: x y */
void writeFile(char s) {
	FILE* file;
	if (s == 'u') {
		fopen_s(&file, "foil_u.txt", "w+");
		if (file == NULL) {
			printf("File open error!\n");
			return;
		}
		for (short i = 0; i <= N; i++)
			fprintf(file, "%lf %lf\n", x_u[i], y_u[i]);
		fclose(file);
	}
	else if (s == 'l') {
		fopen_s(&file, "foil_l.txt", "w+");
		if (file == NULL) {
			printf("File open error!\n");
			return;
		}
		for (short i = 0; i <= N; i++)
			fprintf(file, "%lf %lf\n", x_l[i], y_l[i]);
		fclose(file);
	}
}
int main() {
	generate_foil();
	writeFile('u');
	writeFile('l');
	printf("Write files succeeded!\n");
	draw();
	return 0;
}
