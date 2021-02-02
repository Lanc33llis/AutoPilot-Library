#include <vector>
#include <iostream>
#include <cmath>

const auto PI = 3.141592653589793238462643383279502884L;

void hermite_cubic_to_power_cubic(double x1, double f1, double d1, double x2,
	double f2, double d2, double* c0, double* c1, double* c2, double* c3)

	//****************************************************************************80
	//  Purpose:
	//
	//    HERMITE_CUBIC_TO_POWER_CUBIC converts a Hermite cubic to power form.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    13 February 2011
	//
	//  Author:
	//
	//    John Burkardt
	//
	//  Reference:
	//
	//    Fred Fritsch, Ralph Carlson,
	//    Monotone Piecewise Cubic Interpolation,
	//    SIAM Journal on Numerical Analysis,
	//    Volume 17, Number 2, April 1980, pages 238-246.
	//
	//  Parameters:
	//
	//    Input, double X1, F1, D1, the left endpoint, function value
	//    and derivative.
	//
	//    Input, double X2, F2, D2, the right endpoint, function value
	//    and derivative.
	//
	//    Output, double *C0, *C1, *C2, *C3, the power form of the polynomial.
	//
{
	double df;
	double h;

	h = x2 - x1;
	df = (f2 - f1) / h;
	//
	//  Polynomial in terms of X - X1:
	//
	*c0 = f1;
	*c1 = d1;
	*c2 = -(2.0 * d1 - 3.0 * df + d2) / h;
	*c3 = (d1 - 2.0 * df + d2) / h / h;
	//
	//  Shift polynomial to X.
	//
	*c2 = *c2 - x1 * *c3;
	*c1 = *c1 - x1 * *c2;
	*c0 = *c0 - x1 * *c1;
	*c2 = *c2 - x1 * *c3;
	*c1 = *c1 - x1 * *c2;
	*c2 = *c2 - x1 * *c3;

	return;
}

struct Waypoint
{
    double X, Y, Angle;
    Waypoint(double x, double y, double angleInDegrees) : X(x), Y(y), Angle(angleInDegrees) {}
    Waypoint() : Waypoint(0, 0, 0) {}
};

typedef std::vector<Waypoint> Path;

struct StandardCubicFunction
{
    double A, B, C, D;
    double valueAt(double x) 
    {
        return (A * x * x * x) + (B * x * x) + (C * x)+ D;
    }
    StandardCubicFunction(double a, double b, double c, double d) : A(a), B(b), C(c), D(d) {}
    StandardCubicFunction() : StandardCubicFunction(0, 0, 0, 0) {}
};

struct Spline
{
    Waypoint point1, point2;
    StandardCubicFunction function;
    Spline(Waypoint p1, Waypoint p2, double A, double B, double C, double D) : point1(p1), point2(p2), function(StandardCubicFunction(A, B, C, D)) {}
    Spline(Waypoint p1, Waypoint p2, StandardCubicFunction func) : point1(p1), point2(p2), function(func) {}
};

class Trajectory
{
    std::vector<Spline> curve;
    public:
    Trajectory(std::vector<Spline> splines) : curve(splines) {}
};

double Angle2Deriv(double AngleInDegrees)
{
    double a = AngleInDegrees;
    if (a == 90 || a == 270)
    {
        return NAN;
    }
    else
    {
        return tan(a * (3.141592653589793238462643383279502884197 / 180));
    }
}

Spline HermiteFinder(Waypoint PointOne, Waypoint PointTwo)
{
    // p(x) = c0 + c1 * x + c2 * x^2 + c3 * x^3
    // p(x) = c3 * x^3 + c2 * x^2 + c1 * x + c0
    double *A0 = new double(0), *A1 = new double(0), *A2 = new double(0), *A3 = new double(0);
    hermite_cubic_to_power_cubic(PointOne.X, PointOne.Y, Angle2Deriv(PointOne.Angle), PointTwo.X, PointTwo.Y, Angle2Deriv(PointTwo.Angle), A0, A1, A2, A3);
    return Spline(PointOne, PointTwo, *A3, *A2, *A1, *A0);
}

//uses arc length formula to find distance
double ArcLengthDistance(Spline TheSplineFunction)
{   //deriv = 3ax^2 + 2bx + c
    double Ax = 3 * TheSplineFunction.function.A, Bx = 2 * TheSplineFunction.function.B, C = TheSplineFunction.function.C;
    return (TheSplineFunction.point1.X - TheSplineFunction.point2.X) * (sqrt(1 + pow(pow(Ax, 2) + Bx + C, 2)));
}

//Time given Spline Function and Jerk
double TimeGivenSFJ(Spline TheSplineFunction, double Jerk)
{
    /*

    v(t) = interval [0, t] a(t)dt=1/2Jmax*t^2

    x(t) = interval [0, t] v(t)dt=1/2Jmax*t^3

    thus

    time for end position is:

    t = (6*end position/Jmax)^1/3

    and time for maximum speed

    t = (2*velocitymax/Jmax)^1/2

    */

    double i = ArcLengthDistance(TheSplineFunction);
    return cbrt(6 * ArcLengthDistance(TheSplineFunction) / Jerk);
}

struct Segment
{
    double A, B, C, time;
    Spline spline, xSpline, ySpline;
    double Velocity(double Seconds); 
    double Acceleration(double Seconds); 
    double Jerk(double Seconds);
    Segment(Spline normFunc, Spline xfunc, Spline yfunc, double a, double b, double c, double time) : function(normFunc), xFunction(xfunc), yFunction(yfunc), A(a), B(b), C(c), time(time) {}
    Segment(Spline func, double jerk) : Segment(GenerateSegment(func, jerk)) {}
    private:
        friend Segment GenerateSegment(Spline Function, double Jerk);
};

Segment GenerateSegment(Spline Function, double Jerk)
{
    double B0, B1, B2, Time;
    Time = TimeGivenSFJ(Function, Jerk);
    Spline XFunction = HermiteFinder(Waypoint{ 0, 0, Function.point1.Angle }, Waypoint{ Time, Function.point1.X - Function.point1.X, Function.point2.Angle });
    Spline YFunction = HermiteFinder(Waypoint{ 0, 0, Function.point2.Angle }, Waypoint{ Time, Function.point2.Y - Function.point1.Y, Function.point2.Angle });
    B0 = Function.function.A; B1 = Function.function.B; B2 = Function.function.C;
    return Segment(Function, XFunction, YFunction, B0, B1, B2, Time);
}

