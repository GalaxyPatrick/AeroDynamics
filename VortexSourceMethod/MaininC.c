#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.1415926
#define MAXNUM 50
#define N 8
#define SIGN(x) ((x > 0) - (x < 0))
typedef struct Point Point;
struct Point {
	double x;
	double y;
};
/*采用弧度制*/
const double Velocity_Infinite = 2.0;	//来流速度
const double Radius = 1.0;				//圆柱半径

double tan2cos(double slope);
double tan2sin(double slope);
double func_x(double x, void* Params);
double func_y(double y, void* Params);
Point IntegralSource(Point* PointSet, int n, Point ControlPoint, double* n_value, double* t_value);
void SolveSourceIntensity(double** Coefficient, double* B, double* X);
double* CalculateTanVelocity(double* V_Const, double** Coefficient, double* Source);

int main(void) {
	FILE* file;
	double** n_IntensityCoefficient = (double**)malloc(sizeof(double*) * N);//点源强度系数矩阵
	double** t_IntensityCoefficient = (double**)malloc(sizeof(double*) * N);//切向强度系数矩阵
	double* IntensitySource = (double*)malloc(sizeof(double) * N);		//点源强度
	double* NormalVelocity = (double*)malloc(sizeof(double) * N);		//自由来流法向速度
	double* TangentialVelocity = (double*)malloc(sizeof(double) * N);	//面元来流切向速度
	double* V_t = (double*)malloc(N * sizeof(double));					//面元控制点切向速度

	double InitialAngle = PI / N * 2;		//起始点角度
	Point* CylinderPointSet = (Point*) malloc(N * sizeof(Point));		//面元端点
	Point* CylinderPointControl = (Point*) malloc(N * sizeof(Point));	//面元中点
	Point* tn_PointVelocity = (Point*)malloc(N * sizeof(Point));
	if (CylinderPointControl == NULL || CylinderPointSet == NULL || n_IntensityCoefficient == NULL ||
		NormalVelocity == NULL || t_IntensityCoefficient == NULL) {
		printf("Memory allocation failed!\n");
		return -1;
	}
	else {
		for (int i = 0; i < N; i++) {
			*(n_IntensityCoefficient + i) = (double*)malloc(sizeof(double) * N);
			*(t_IntensityCoefficient + i) = (double*)malloc(sizeof(double) * N);
		}
	}
	/*计算面元控制点坐标*/
	for (short i = 0; i < N; i++){
		(*(CylinderPointSet + i)).x = cos(PI / N + i * 2 * PI / N);
		(*(CylinderPointSet + i)).y = sin(PI / N + i * 2 * PI / N);
		(*(CylinderPointControl + i)).x = cos(InitialAngle + i * 2 * PI / N);
		(*(CylinderPointControl + i)).y = sin(InitialAngle + i * 2 * PI / N);
	}

	for (short i = 0; i < N; i++) {
		*(tn_PointVelocity + i) = IntegralSource(CylinderPointSet, i, CylinderPointControl[i],
			n_IntensityCoefficient[i], t_IntensityCoefficient[i]);
		NormalVelocity[i] = - Velocity_Infinite * (*(tn_PointVelocity + i)).x;
		TangentialVelocity[i] = Velocity_Infinite * (*(tn_PointVelocity + i)).y;
		//printf("%lf %lf\n", (*(tn_PointVelocity + i)).x, (*(tn_PointVelocity + i)).y);
	}
	SolveSourceIntensity(n_IntensityCoefficient, NormalVelocity, IntensitySource);
	V_t = CalculateTanVelocity(TangentialVelocity, t_IntensityCoefficient, IntensitySource);

	/*输出,写文件*/
	printf("Begin to write file.\n");
	fopen_s(&file, "result.txt", "w+");
	if (file == NULL) {
		printf("File open error!\n");
		exit(0);
	}
	for (short i = 0; i < N; i++) {
		//printf("Cp%d = %lf\n", i, 1 - pow(V_t[i] / Velocity_Infinite, 2));
		double tmp_x = atan2(CylinderPointControl[i].y, CylinderPointControl[i].x);
		if (tmp_x < 0) tmp_x = 2 * PI + tmp_x;
		fprintf(file, "%lf %lf\n", tmp_x, 1 - pow(V_t[i] / Velocity_Infinite, 2));
	}
	fclose(file);
	printf("File written succeeded.\n");
	free(CylinderPointSet);
	free(CylinderPointControl);
	return 0;
}
double* CalculateTanVelocity(double* V_Const, double** Coefficient, double* Source) {
	double TempCoefficient[MAXNUM * MAXNUM];
	gsl_matrix_view A;
	gsl_vector_view X;
	gsl_vector_view V;
	double* velocity_result = (double*)malloc(N * sizeof(double));

	for (short i = 0; i < N * N; i++) {
		TempCoefficient[i] = Coefficient[i / N][i % N];
	}
	A = gsl_matrix_view_array(TempCoefficient, N, N);	//转换为GSL内部矩阵类型
	X = gsl_vector_view_array(Source, N);
	V = gsl_vector_view_array(V_Const, N);
	gsl_blas_dgemv(CblasNoTrans, 1.0, &A, &X, 1.0, &V);
	printf("v_t = \n");
	gsl_vector_fprintf(stdout, &V, "%g");
	for (short i = 0; i < N; i++) {
		velocity_result[i] = V.vector.data[i];
	}
	return velocity_result;
}
/*求解点源强度矩阵*/
void SolveSourceIntensity(double** Coefficient, double* B, double* X) {
	double TempCoefficient[MAXNUM * MAXNUM];
	gsl_matrix_view A;
	gsl_vector_view b = gsl_vector_view_array(B, N);
	gsl_vector* x = gsl_vector_alloc(N);
	int s;
	//printf("tempCoefficient = \n");
	for (short i = 0; i < N * N; i++) {
		TempCoefficient[i] = Coefficient[i / N][i % N];
		//printf("%lf ", TempCoefficient[i]);
		//if ((i + 1) % N == 0) printf("\n");
	}
	A = gsl_matrix_view_array(TempCoefficient, N, N);	//转换为GSL内部矩阵类型
	gsl_permutation* p = gsl_permutation_alloc(N);
	gsl_linalg_LU_decomp(&A.matrix, p, &s);
	gsl_linalg_LU_solve(&A.matrix, p, &b.vector, x);
	
	printf("x = \n");
	gsl_vector_fprintf(stdout, x, "%g");
	gsl_permutation_free(p);
	for (int i = 0; i < N; i++) {
		X[i] = x->data[i];
	}
	printf("SourceIntensity calculated succeeded!\n");
	gsl_vector_free(x);
}

/*PointSet提供场源面点集，n表示积分场点序号，ControlPoint提供场源点
n_value返回法向场源强度方程系数行向量
IntegralSource以单位向量形式返回第n个面的法向量（余弦、正弦）*/
Point IntegralSource(Point* PointSet, int n, Point ControlPoint, double* n_value, double* t_value) {
	/*积分（含奇点）迭代空间*/
	//gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
	/*被积函数结构体*/
	gsl_function func_x1, func_y1;
	/*积分参数列表*/
	double IntegralParams[4] = {ControlPoint.x, ControlPoint.y, 0, 0};
	Point* a;
	Point* b;
	double _normal_velocity_coefficient;	//法向速度系数
	double _tan_velocity_coefficient;		//切向速度系数
	double res, partialX, partialY;
	size_t neval;
	double slope;					//场源面斜率
	double slope_surface;			//场点面斜率
	double cos_n, sin_n;			//场点面法向量
	double cos_t, sin_t;			//场点面方向向量
	Point PointVector;				
	/*----计算场点面法向量、方向向量----*/
	a = PointSet + n;
	b = PointSet + (n + 1) % N;
	slope_surface = (a->y - b->y) / (a->x - b->x);
	cos_t = SIGN(b->y - a->y) * tan2cos(slope_surface);
	sin_t = SIGN(b->x - a->x) * tan2sin(slope_surface);
	cos_n = tan2cos(-1 / slope_surface);
	sin_n = tan2sin(-1 / slope_surface);
	if (((b->x - a->x) * sin_n - (b->y - a->y) * cos_n) > 0) {
		cos_n = -cos_n;
		sin_n = -sin_n;
	}
	PointVector.x = cos_n;
	PointVector.y = sin_n;
	//遍历场源面
	for (short i = 0; i < N; i++) {
		if (i == n) {				//积分场点在场源上
			_normal_velocity_coefficient = 0.5;		//奇点
			_tan_velocity_coefficient = 0.0;
		}
		else {
			/*----场源面方程y = slope x + b----*/
			a = PointSet + i;
			b = PointSet + (i + 1) % N;
			slope = (a->y - b->y) / (a->x - b->x);
			double c = sqrt(1 + pow(slope, 2)) / 4 / PI;
			double d = sqrt(1 + pow(1 / slope, 2)) / 4 / PI;
			/*----IntegralParams[]----
			[0][1]: 场点坐标
			[2]: 场源面斜率
			[3]: 场源面截距*/
			/*----x方向导数----*/
			IntegralParams[2] = slope;
			IntegralParams[3] = a->y - slope * a->x;
			func_x1.function = &func_x;
			func_x1.params = IntegralParams;
			gsl_integration_qng(&func_x1, a->x, b->x, 0, 1e-7,
				&partialX, &res, &neval);

			/*----y方向导数----*/
			IntegralParams[2] = 1 / slope;
			IntegralParams[3] = -IntegralParams[3] / slope;
			func_y1.function = &func_y;
			func_y1.params = IntegralParams;
			gsl_integration_qng(&func_y1, a->y, b->y, 0, 1e-7,
				&partialY, &res, &neval);

			_normal_velocity_coefficient = c * partialX * cos_n + d * partialY * sin_n;
			_tan_velocity_coefficient = c * partialX * cos_t + d * partialY * sin_t;
		}
		*(n_value + i) = _normal_velocity_coefficient;			//场源强度系数矩阵(第n行向量)
		*(t_value + i) = _tan_velocity_coefficient;
	}
	//返回面法向量余弦
	//return cos_n;
	return PointVector;
}
double tan2cos(double slope)
{
	return 1 / sqrt(1 + slope * slope);
}
double tan2sin(double slope)
{
	return slope / sqrt(1 + slope * slope);
}
double func_x(double x, void* Params) {
	Point dPoint;
	/*----参数列表解引用----*/
	double k = ((double*)Params)[2];
	double b = ((double*)Params)[3];
	dPoint.x = ((double*)Params)[0];
	dPoint.y = ((double*)Params)[1];

	double r2 = pow(x - dPoint.x, 2) + pow(k * x + b - dPoint.y, 2);
	double func_x = 2 * (x - dPoint.x) / r2;
	return func_x;
}
double func_y(double y, void* Params) {
	Point dPoint;
	/*----参数列表解引用----*/
	double k = ((double*)Params)[2];
	double b = ((double*)Params)[3];
	dPoint.x = ((double*)Params)[0];
	dPoint.y = ((double*)Params)[1];

	double r2 = pow(k * y + b - dPoint.x, 2) + pow(y - dPoint.y, 2);
	double func_y = 2 * (y - dPoint.y) / r2;
	return func_y;
}