#include <iostream>
using namespace std;
const double Velocity_Infinite = 2.0;
class Point {
public:
	Point(double a = 0, double b = 0):x(a), y(b){}
	double getX();
	double getY();
	friend ostream& operator << (ostream& output, Point& p);
protected:
	double x;
	double y;
};

double Point::getX()   
{
	return x;
}
double Point::getY()
{
	return y;
}
ostream& operator << (ostream& output, Point& p)
{
	output << "(" << p.x << "," << p.y << ">" << endl;
	return output;
}
int main() {
	Point p1(2.5, 3.7);
	std::cout << "hello" << endl;
	return 0;
}