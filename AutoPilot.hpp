#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>

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

double truncate(double input, size_t accuracy)
{
    bool isNeg = false;
    if (input < 0) 
    {
        isNeg = true;
    }
    std::string a = std::to_string(input);
    std::string c = std::string();
    if (isNeg) {
        c = a.substr(0, accuracy + 3);
    }
    else {
        c = a.substr(0, accuracy + 2);
    }
    return stod(c);
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
        return (A * x * x * x) + (B * x * x) + (C * x) + D;
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
    Spline() : Spline(HermiteFinder(Waypoint(0, 0, 0), Waypoint(0, 0, 0))) {}
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
    double* A0 = new double(0), * A1 = new double(0), * A2 = new double(0), * A3 = new double(0);
    hermite_cubic_to_power_cubic(PointOne.X, PointOne.Y, Angle2Deriv(PointOne.Angle), PointTwo.X, PointTwo.Y, Angle2Deriv(PointTwo.Angle), A0, A1, A2, A3);
    return Spline(PointOne, PointTwo, *A3, *A2, *A1, *A0);
}

//Spline HermiteFinder(Waypoint PointOne, Waypoint PointTwo)
//{
//    // p(x) = c0 + c1 * x + c2 * x^2 + c3 * x^3
//    // p(x) = c3 * x^3 + c2 * x^2 + c1 * x + c0
//    double* A0 = new double(0), * A1 = new double(0), * A2 = new double(0), * A3 = new double(0);
//    hermite_cubic_to_power_cubic(PointOne.X, PointOne.Y, Angle2Deriv(PointOne.Angle), PointTwo.X, PointTwo.Y, Angle2Deriv(PointTwo.Angle), A0, A1, A2, A3);
//    return Spline(PointOne, PointTwo, *A3, *A2, *A1, *A0);
//    delete A0;
//    delete A1;
//    delete A2;
//    delete A3;
//}

//uses arc length formula to find distance
double ArcLengthDistance(Spline Function, size_t Accuracy = 30)
{   //deriv = 3ax^2 + 2bx + c
    //http://www2.cs.uregina.ca/~anima/408/Notes/MotionControl/AnalyticApproach.htm
    /*
    A = 9(Ax^2+Ay^2)
    B = 12(Ax*Bx+Ay*By)
    C = 6(Ax*Cx+Ay*Cy)+4*(Bx^2+By^2)
    D = 4(Bx*Cx+By*Cy)
    E = Cx^2+Cy^2
    */
    double x1 = Function.point1.X;
    double x2 = Function.point2.X;

    Spline XSpline = HermiteFinder(Waypoint( 0, Function.point1.X, 45), Waypoint( Function.point2.X - Function.point1.X, Function.point2.X, 45));
    Spline YSpline = HermiteFinder(Waypoint( 0, Function.point1.Y, Function.point1.Angle ), Waypoint( Function.point2.X - Function.point1.X, Function.point2.Y, Function.point2.Angle ));
    auto f = [XSpline, YSpline](double u)
    {
        auto xfunc = XSpline.function;
        auto yfunc = YSpline.function;
        double Ax = xfunc.A;
        double Bx = xfunc.B;
        double Cx = xfunc.C;
        double Dx = xfunc.D;
        double Ay = yfunc.A;
        double By = yfunc.B;
        double Cy = yfunc.C;
        double Dy = yfunc.D;
        double A = 9*(Ax*Ax+Ay*Ay);
        double B = 12*(Ax*Bx+Ay*By);
        double C = 6*(Ax*Cx+Ay*Cy)+4*(Bx*Bx+By*By);
        double D = 4*(Bx*Cx+By*Cy);
        double E = Cx*Cx+Cy*Cy;

        return sqrt(A*pow(u, 4) + B*pow(u, 3) + C*pow(u, 2) + D*u + E);
    };

    double arcLength = 0;
    double h = (x2 - x1) / Accuracy;
    for (size_t i = 1; i <= Accuracy; i++)
    {
        arcLength += h/2 * (f(x1+(i-1)*h) + f(x1+i*h));
    }
    return abs(arcLength);
}

double ArcLengthDistance(Spline theSplineFunction, double lowerLimit, double upperLimit)
{
    theSplineFunction.point2.X = upperLimit;
    theSplineFunction.point1.X = lowerLimit;
    return ArcLengthDistance(theSplineFunction);
}

// double DistanceFromPointUsingArcLength(Spline theSpline, double x, double distance)
// {
//     //distance = (upper - lower) * √(1 + f'(x)²)
//     double Ax = 3 * theSpline.function.A, Bx = 2 * theSpline.function.B, C = theSpline.function.C;
//     //distance = (theSpline.point2.X - theSpline.point1.X) * √(1 + (Ax*Ax + Bx + c)²)
// }

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
    //TIME STARTS AT 0 THAT'S WHY FUNCTION STARTS at 0, 0
    double A, B, C, D, Time;
    Time = TimeGivenSFJ(Function, Jerk);
    Spline XFunction = HermiteFinder(Waypoint( 0, 0, Function.point1.Angle ), Waypoint( Time, Function.point2.X - Function.point1.X, Function.point2.Angle ));
    Spline YFunction = HermiteFinder(Waypoint( 0, 0, Function.point1.Angle ), Waypoint( Time, Function.point2.Y - Function.point1.Y, Function.point2.Angle ));
    A = Function.function.A; B = Function.function.B; C = Function.function.C; D = Function.function.D;
    return Segment(Function, XFunction, YFunction, A, B, C, D, Time);
}

//these values determine how many segments are created.
//USING_DISTANCE_PARTION determines if a new segment is created based off distance or a equal division
//TANK_CONFIG_PARTION_VALUE represents the distance or how many divisions to create
bool USING_DISTANCE_PARTION = true;
size_t TANK_CONFIG_PARTION_VALUE = .25;

typedef std::vector<Segment> Trajectory;
class TankConfig
{
    Curve curve;

    Trajectory leftTrajectory;
    Trajectory rightTrajectory;

    friend void createDesmosGraph(TankConfig config, std::string fileName, std::string htmlFileLocation);

    public:
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
            std::cout << "Right Segment " << i << "'s Values: " << rightTrajectory[i].A << " " << rightTrajectory[i].B << " " << rightTrajectory[i].C << " " << rightTrajectory[i].D << "\n";
        }
    }

    TankConfig(Curve curve, double widthBetweenWheels, double jerk)
    {
        TankConfig::curve = curve;
        leftTrajectory = Trajectory();
        rightTrajectory = Trajectory();
        double centerToWheels = widthBetweenWheels / 2;
        for (Spline s : curve)
        {
            auto x1 = s.point1.X;
            auto y1 = s.point1.Y;
            auto slope1 = truncate((3 * s.function.A * x1 * x1) + (2 * s.function.B * x1) + (s.function.C), 10);
            auto normal1 = -1 / slope1;
            auto x2 = s.point2.X;
            auto y2 = s.point2.Y;
            auto slope2 = truncate((3 * s.function.A * x2 * x2) + (2 * s.function.B * x2) + (s.function.C), 10);
            auto normal2 = -1 / slope2;

            auto pointSlopeForm = [](double x1, double y1, double m)
            {
                return [x1, y1, m](double x) -> double {
                    return (m * (x - x1) + y1);
                };
            };

            auto normalLine1 = pointSlopeForm(x1, y1, normal1);
            auto normalLine2 = pointSlopeForm(x2, y2, normal2);

            double nx1, nx2, nx3, nx4, ny1, ny2, ny3, ny4;

            if (normal1 == NAN || normal1 == -INFINITY || normal1 == -NAN || normal1 == INFINITY)
            {
                nx1 = x1;
                nx3 = x1;
                if (sin(s.point1.Angle * DEGREES2RADIANS) > 0)
                {
                    ny1 = y1 + centerToWheels;
                    ny3 = y1 - centerToWheels;
                }
                else
                {
                    ny1 = y1 - centerToWheels;
                    ny3 = y1 + centerToWheels;
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
            if (normal2 == NAN || normal2 == -INFINITY || normal2 == -NAN || normal2 == INFINITY)
            {
                nx2 = x2;
                nx4 = x2;
                if (sin(s.point2.Angle * DEGREES2RADIANS) > 0)
                {
                    ny2 = y2 + centerToWheels;
                    ny4 = y2 - centerToWheels;
                }
                else
                {
                    ny2 = y2 - centerToWheels;
                    ny4 = y2 + centerToWheels;
                }
            }
            else
            {
                if (sin(s.point2.Angle * DEGREES2RADIANS) > 0)
                {
                    nx2 = x2 + (centerToWheels / sqrt(1 + (normal2 * normal2)));
                    nx4 = x2 - (centerToWheels / sqrt(1 + (normal2 * normal2)));
                }
                else
                {
                    nx2 = x2 - (centerToWheels / sqrt(1 + (normal2 * normal2)));
                    nx4 = x2 + (centerToWheels / sqrt(1 + (normal2 * normal2)));
                }
                ny2 = normalLine2(nx2);
                ny4 = normalLine2(nx4);
            }

            
            Waypoint np1 = Waypoint(nx1, ny1, s.point1.Angle), np2 = Waypoint(nx2, ny2, s.point2.Angle),
            np3 = Waypoint(nx3, ny3, s.point1.Angle), np4 = Waypoint(nx4, ny4, s.point2.Angle);

            auto s1 = Spline(np1, np2);
            auto s2 = Spline(np3, np4);

            leftTrajectory.push_back(Segment(s1, jerk));
            rightTrajectory.push_back(Segment(s2, jerk));
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

//creates html file that displays desmos graph, relative location files only work if your workspace location is correct
//htmlFileLocation overrides the file name
void createDesmosGraph(TankConfig config, std::string fileName = "graph.html", std::string htmlFileLocation = "")
{
    using namespace std::string_literals;
    using namespace std;
    ofstream file;

    if (htmlFileLocation.empty())
    {
        file.open("./" + fileName);
    }
    else
    {
        file.open(htmlFileLocation);
    }
    file << "<!DOCTYPE html><html><head></head><style>html, body{margin: 0px; width: 100%; height: 100%}</style><body>";
    file << u8R"(<script src="https://www.desmos.com/api/v1.5/calculator.js?apiKey=dcb31709b452b1cf9dc26972add0fda6"></script>)"s;


    file << u8R"(<div id="calculator" style="width: 100%; height: 100%;"></div>)"s;

    file << "<script>\n";
    file << u8R"(var elt = document.getElementById('calculator');)"s << "\n" << "var calculator = Desmos.GraphingCalculator(elt);";

    size_t precision = 10;

    for (size_t i = 0; i < config.leftTrajectory.size(); i++)
    {
        auto s = config.leftTrajectory[i];
        file << u8R"(calculator.setExpression({id: 'leftGraph)"s;
        file << i;
        file << u8R"(', color: Desmos.Colors.RED, latex: 'y=)"s;
        file << fixed << setprecision(precision) << s.A << "x^3 + " << setprecision(precision) << fixed << s.B << "x^2 + " << setprecision(precision) << fixed << s.C << "x + " << setprecision(precision) << fixed << s.D;
        file << "  \\\\left\\\\{" << setprecision(6) << min(s.spline.point1.X, s.spline.point2.X) << "\\\\le x \\\\le" << setprecision(6) << max(s.spline.point1.X, s.spline.point2.X) << "\\\\right\\\\}";
        file << u8R"('});)"s;
    }

    for (size_t i = 0; i < config.rightTrajectory.size(); i++)
    {
        auto s = config.rightTrajectory[i];
        file << u8R"(calculator.setExpression({id: 'rightGraph)"s;
        file << i;
        file << u8R"(', color: Desmos.Colors.BLUE, latex: 'y=)"s;
        file << setprecision(precision) << fixed << s.A << "x^3 + " << setprecision(precision) << fixed << s.B << "x^2 + " << setprecision(precision) << fixed << s.C << "x + " << setprecision(precision) << fixed << s.D;
        file << "  \\\\left\\\\{" << setprecision(6) << min(s.spline.point1.X, s.spline.point2.X) << "\\\\le x \\\\le" << setprecision(6) << max(s.spline.point1.X, s.spline.point2.X) << "\\\\right\\\\}";
        file << u8R"('});)"s;
    }

    //cout << "curve size: " << config.curve.size() << "\n";
    for (size_t i = 0; i < config.curve.size(); i++)
    {
        file << u8R"(calculator.setExpression({id: 'test)"s << i << "'";
        file << u8R"(, color: Desmos.Colors.BLACK, latex: 'y=)"s;
        file << setprecision(precision) << fixed << config.curve[i].function.A << "x^3 + " << setprecision(precision) << fixed << config.curve[i].function.B << "x^2 + " << setprecision(precision) << fixed << config.curve[i].function.C << "x + " << setprecision(precision) << fixed << config.curve[i].function.D;
        file << "  \\\\left\\\\{" << setprecision(6) << min(config.curve[i].point1.X, config.curve[i].point2.X) << "\\\\le x \\\\le" << setprecision(6) << max(config.curve[i].point1.X, config.curve[i].point2.X) << "\\\\right\\\\}";
        file << u8R"('});)"s;
    }
    file << "</script>";

    file.close();
}

void createDesmosGraph(Curve curve, std::string fileName = "graph.html", std::string htmlFileLocation = "")
{
    using namespace std::string_literals;
    using namespace std;
    ofstream file;

    if (htmlFileLocation.empty())
    {
        file.open("./" + fileName);
    }
    else
    {
        file.open(htmlFileLocation);
    }
    file << "<!DOCTYPE html><html><head></head><style>html, body{margin: 0px; width: 100%; height: 100%}</style><body>";
    file << u8R"(<script src="https://www.desmos.com/api/v1.5/calculator.js?apiKey=dcb31709b452b1cf9dc26972add0fda6"></script>)"s;


    file << u8R"(<div id="calculator" style="width: 100%; height: 100%;"></div>)"s;

    file << "<script>\n";
    file << u8R"(var elt = document.getElementById('calculator');)"s << "\n" << "var calculator = Desmos.GraphingCalculator(elt);";

    size_t precision = 10;

    for (size_t i = 0; i < curve.size(); i++)
    {
        file << u8R"(calculator.setExpression({id: 'test)"s << i << "'";
        file << u8R"(, color: Desmos.Colors.BLACK, latex: 'y=)"s;
        file << setprecision(precision) << fixed << curve[i].function.A << "x^3 + " << setprecision(precision) << fixed << curve[i].function.B << "x^2 + " << setprecision(precision) << fixed << curve[i].function.C << "x + " << setprecision(precision) << fixed << curve[i].function.D;
        file << "  \\\\left\\\\{" << setprecision(6) << min(curve[i].point1.X, curve[i].point2.X) << "\\\\le  x \\\\le" << setprecision(6) << max(curve[i].point1.X, curve[i].point2.X) << "\\\\right\\\\} ";
        file << u8R"('});)"s;
    }
    file << "</script>";

    file.close();
}