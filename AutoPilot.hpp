#include <vector>
#include <iostream>
#include <cmath>

const auto PI = 3.141592653589793238462643383279502884L;
const auto DEGREES2RADIANS = (PI / 180);
const auto RADIANS2DEGREES = (180 / PI);

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

struct Spline;
Spline HermiteFinder(Waypoint PointOne, Waypoint PointTwo);
struct Spline
{
    Waypoint point1, point2;
    StandardCubicFunction function;
    Spline(Waypoint p1, Waypoint p2, double A, double B, double C, double D) : point1(p1), point2(p2), function(StandardCubicFunction(A, B, C, D)) {}
    Spline(Waypoint p1, Waypoint p2, StandardCubicFunction func) : point1(p1), point2(p2), function(func) {}
    Spline(Waypoint p1, Waypoint p2) : Spline(HermiteFinder(p1, p2)) {}
    Spline() : Spline(HermiteFinder(Waypoint(0,0,0), Waypoint(0,0,0))) {}
};

typedef std::vector<Spline> Curve;

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
    return abs(TheSplineFunction.point2.X - TheSplineFunction.point1.X) * (sqrt(1 + pow(pow(Ax, 2) + Bx + C, 2)));
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

    return cbrt(6 * ArcLengthDistance(TheSplineFunction) / Jerk);
}

struct Segment;
Segment GenerateSegment(Spline Function, double Jerk);
struct Segment
{
    double A, B, C, D, time;
    Spline spline, xSpline, ySpline;
    double Velocity(double Seconds) 
    {
		if (Seconds > time)
		{
			return NAN;
		}
		else
		{
			double Xa = xSpline.function.A, Xb = xSpline.function.B, Xc = xSpline.function.C;
			double Ya = ySpline.function.A, Yb = ySpline.function.B, Yc = ySpline.function.C;
			return sqrt(pow((3 * Xa * pow(Seconds, 2)) + (2 * Xb * Seconds) + Xc, 2) + pow((3 * Ya * pow(Seconds, 2)) + (2 * Yb * Seconds) + Yc, 2));
		}
    }

    double Acceleration(double Seconds) 
    {
        if (Seconds > time)
		{
			return NAN;
		}
		else
		{
			double Xa = xSpline.function.A, Xb = xSpline.function.B;
			double Ya = ySpline.function.A, Yb = ySpline.function.B;
			// return sqrt(pow((3 * Xa * pow(Seconds, 2)) + (2 * Xb * Seconds) + Xc, 2) + pow((3 * Ya * pow(Seconds, 2)) + (2 * Yb * Seconds) + Yc, 2));
            return sqrt(pow(((6 * Xa * Seconds) + (2 * Xb)), 2) + pow(((6 * Ya * Seconds) + (2 * Yb)), 2));
		}
    }

    double Jerk(double Seconds) 
    {
        if (Seconds > time)
		{
			return NAN;
		}
		else
		{
			double Xa = xSpline.function.A;
			double Ya = ySpline.function.A;
			// return sqrt(pow((3 * Xa * pow(Seconds, 2)) + (2 * Xb * Seconds) + Xc, 2) + pow((3 * Ya * pow(Seconds, 2)) + (2 * Yb * Seconds) + Yc, 2));
            return sqrt(pow(6 * Xa, 2) + pow(6 * Ya, 2));
		}
    }

    Segment(Spline normFunc, Spline xfunc, Spline yfunc, double a, double b, double c, double d, double time) : spline(normFunc), xSpline(xfunc), ySpline(yfunc), A(a), B(b), C(c), D(d), time(time) {}
    Segment(Spline func, double jerk) : Segment(GenerateSegment(func, jerk)) {}
    Segment() : A(0), B(0), C(0), D(0), time(0), spline(Spline()), xSpline(Spline()), ySpline(Spline()) {}
};

Segment GenerateSegment(Spline Function, double Jerk)
{
    double A, B, C, D, Time;
    Time = TimeGivenSFJ(Function, Jerk);
    Spline XFunction = HermiteFinder(Waypoint{ 0, 0, Function.point1.Angle }, Waypoint{ Time, Function.point1.X - Function.point1.X, Function.point2.Angle });
    Spline YFunction = HermiteFinder(Waypoint{ 0, 0, Function.point2.Angle }, Waypoint{ Time, Function.point2.Y - Function.point1.Y, Function.point2.Angle });
    A = Function.function.A; B = Function.function.B; C = Function.function.C; D = Function.function.D;
    return Segment(Function, XFunction, YFunction, A, B, C, D, Time);
}

typedef std::vector<Segment> Trajectory;
class TankConfig 
{
    public:
    Trajectory leftTrajectory;
    Trajectory rightTrajectory;

    template<typename SpeedControllerGroup, typename Encoder>
    void run(SpeedControllerGroup leftSide, SpeedControllerGroup rightSide, double jerk, double Kp)
    {
        //u = Κₚe

    }

    void testTrajectory()
    {
        for (size_t i = 0; i < leftTrajectory.size(); i++)
        {
            std::cout << "Left Segment " << i << "'s Values: " << leftTrajectory[i].A << " " << leftTrajectory[i].B << " " << leftTrajectory[i].C << " " << leftTrajectory[i].D << " " << leftTrajectory[i].spline.point1.X << " " << leftTrajectory[i].spline.point1.Y << "\n";        
        }
        for (size_t i = 0; i < rightTrajectory.size(); i++)
        {
            std::cout << "Right Segment " << i << "'s Values: " << rightTrajectory[i].A << " " << rightTrajectory[i].B << " " << rightTrajectory[i].C << " " << rightTrajectory[i].D <<"\n";        
        }
    }

    TankConfig(Curve curve, double widthBetweenWheels, double jerk) 
    {
        leftTrajectory = Trajectory();
        rightTrajectory = Trajectory();
        double centerToWheels = widthBetweenWheels / 2;
        for (Spline s : curve) 
        {
            auto x1 = s.point1.X;
            auto y1 = s.point1.Y;
            auto slope1 = (3 * s.function.A * x1 * x1) + (2 * s.function.B * x1) + (s.function.C);
            auto normal1 = -1 / slope1;
            auto x2 = s.point2.X;
            auto y2 = s.point2.Y;
            auto slope2 =(3 * s.function.A * x2 * x2) + (2 * s.function.B * x2) + (s.function.C);
            auto normal2 = -1 / slope2;

            auto pointSlopeForm = [](double x1, double y1, double m) 
            {
                return [x1, y1, m](double x) -> double {
                    return (m * (x - x1) + y1);
                };
            };
            
            auto normalLine1 = pointSlopeForm(x1, y1, normal1);
            auto normalLine2 = pointSlopeForm(x2, x2, normal2);

            double nx1, nx2, nx3, nx4, ny1, ny2, ny3, ny4;
            
            if (normal1 == NAN || normal1 == -INFINITY) 
            {
                nx1 = x1; 
                nx3 = x1;
                if (sin(s.point1.Angle * DEGREES2RADIANS) > 0)
                {
                    ny1 = y1 - centerToWheels;
                    ny3 = y1 + centerToWheels;
                }
                else 
                {
                    ny1 = y1 + centerToWheels;
                    ny3 = y1 - centerToWheels;
                }
            } 
            else {
                if (sin(s.point1.Angle * DEGREES2RADIANS) > 0)
                {
                    nx1 = x1 + (centerToWheels / sqrt(1 + (normal1 * normal1)));
                    nx3 = x1 - (centerToWheels / sqrt(1 + (normal1 * normal1)));
                }
                else 
                {
                    nx1 = x1 - (centerToWheels / sqrt(1 + (normal1 * normal1)));
                    nx3 = x1 + (centerToWheels / sqrt(1 + (normal1 * normal1)));
                }
                ny1 = normalLine1(nx1);
                ny3 = normalLine1(nx3);
            }
            if (normal2 == NAN || normal2 == -INFINITY)
            {
                nx2 = x2;
                nx4 = x2;
                if (sin(s.point2.Angle * DEGREES2RADIANS) > 0)
                {
                    ny2 = y2 - centerToWheels;
                    ny4 = y2 + centerToWheels;
                }
                else 
                {
                    ny2 = y2 + centerToWheels;
                    ny4 = y2 - centerToWheels;
                }
            } 
            else 
            {
                if (sin(s.point2.Angle * DEGREES2RADIANS) > 0) 
                {
                    nx4 = x2 - (centerToWheels / sqrt(1 + (normal2 * normal2)));
                    nx2 = x2 + (centerToWheels / sqrt(1 + (normal2 * normal2)));
                } 
                else
                {
                    nx4 = x2 + (centerToWheels / sqrt(1 + (normal2 * normal2)));
                    nx2 = x2 - (centerToWheels / sqrt(1 + (normal2 * normal2)));
                }
                ny2 = normalLine2(nx2);
                ny4 = normalLine2(nx4);
            }


            Waypoint np1 = Waypoint(nx1, ny1, s.point1.Angle), np2 = Waypoint(nx2, ny2, s.point2.Angle),
            np3 = Waypoint(nx3, ny3, s.point1.Angle), np4 = Waypoint(nx4, ny4, s.point2.Angle);
            std::cout << "Wp1 " << x1 << " " <<  y1 << " " << s.point1.Angle << "\n";
            std::cout << "Wp2 " << x2 << " " <<  y2 << " " << s.point2.Angle << "\n";
            std::cout << "np1 " << np1.X << " " <<  np1.Y << " " << np1.Angle << "\n";
            std::cout << "np2 " << np2.X << " " <<  np2.Y << " " << np2.Angle << "\n";
            std::cout << "np3 " << np3.X << " " <<  np3.Y << " " << np3.Angle << "\n";
            std::cout << "np4 " << np4.X << " " <<  np4.Y << " " << np4.Angle << "\n";
            // leftSideCurve.push_back(Spline(Waypoint(nx1, ny2, s.point1.Angle), Waypoint(nx2, ny2, s.point1.Angle)));
            // rightSideCurve.push_back(Spline(Waypoint(nx3, ny3, s.point1.Angle), Waypoint(nx4, ny4, s.point2.Angle)));
            leftTrajectory.push_back(Segment(Spline(np1, np2), jerk));
            rightTrajectory.push_back(Segment(Spline(np3, np4), jerk));
        }
    }
};

typedef std::vector<Waypoint> Path;

Curve curveGenerator(Path path) {
    Curve ReturnSpline;
    size_t NumberOfFunctions = path.size() - 1;
    for (int i = 0; i < NumberOfFunctions; i++)
    {
        auto Temp = HermiteFinder(path[i], path[i + 1]);
        ReturnSpline.push_back(Temp);
    }
    return ReturnSpline;
}