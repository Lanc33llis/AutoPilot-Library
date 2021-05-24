#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <initializer_list>

const auto PI = 3.141592653589793238462643383279502884L;
const auto DEGREES2RADIANS = (PI / 180);
const auto RADIANS2DEGREES = (180 / PI);

class Mat
{
    size_t rows, cols;
    std::vector<std::vector<double>> arr;
    std::vector<std::vector<double>> powerArr;
    std::vector<std::vector<double>> varArr;
public:
    Mat(size_t rows, size_t cols, double fill, double pow) : rows(rows), cols(cols)
    {
        arr = std::vector<std::vector<double>>(cols, std::vector<double>(rows, fill));
        powerArr = std::vector<std::vector<double>>(cols, std::vector<double>(rows, pow));
        varArr = std::vector<std::vector<double>>(cols, std::vector<double>(rows, 0));
    }
    Mat(size_t rows, size_t cols, double fill) : Mat(rows, cols, fill, 1) {}
    Mat(size_t rows, size_t cols) : Mat(rows, cols, 0, 1) {}
    Mat(std::initializer_list<std::vector<double>> list)
    {
        arr = list;
        cols = arr.size();
        rows = arr[0].size();

        powerArr = std::vector<std::vector<double>>(cols, std::vector<double>(rows, 1));
        varArr = std::vector<std::vector<double>>(cols, std::vector<double>(rows, 0));
    }
    Mat Mult(Mat b)
    {
        Mat A = *this;
        Mat B = b;

        if (A.cols == B.rows)
        {
            Mat n = Mat(A.rows, B.cols);

            for (size_t i = 0; i < n.rows; i++)
            {
                for (size_t j = 0; j < n.cols; j++)
                {
                    for (size_t k = 0; k < B.rows; k++)
                    {
                        n[j][i] += A[k][i] * B[j][k];
                    }
                }
            }

            return n;
        }
        else
        {
            std::string s = std::string("Invalid Matrix Dimensions, A has ") + std::to_string(A.cols) + std::string(" cols ") + std::string("but B has ") + std::to_string(B.rows) + std::string(" rows");
            std::cerr << s << "\n";
            throw std::logic_error(s);
        }
    }

    void setPowers(Mat a)
    {
        powerArr = a.arr;
    }

    void setVariables(Mat a)
    {
        varArr = a.arr;
    }

    Mat input(double t)
    {
        Mat b(rows, cols);
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < cols; j++)
            {
                if (varArr[j][i] == 1)
                {
                    b[j][i] = std::pow(t, powerArr[j][i]) * arr[j][i];
                }
                else
                {
                    b[j][i] = arr[j][i];
                }
            }
        }
        return b;
    }

    std::vector<double>& operator [](const size_t i)
    {
        return arr[i];
    }

    operator std::vector < std::vector<double> >()
    {
        return arr;
    }

    void print()
    {
        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < cols; j++)
            {
                std::cout << arr[j][i] << " ";
            }
            std::cout << "\n";
        }
    }

    static Mat InverseBeizerBasis()
    {
        //Mat t(4, 4);
        //t[0][0] = 1;
        //t[0][1] = 1;
        //t[0][2] = 1;
        //t[0][3] = 1;

        //t[1][1] = 1 / 3;
        //t[1][2] = 2 / 3;
        //t[1][3] = 1;

        //t[2][2] = 1 / 3;
        //t[2][3] = 1;

        //t[3][3] = 1;

        Mat t = { {1, 1, 1, 1}, {0, (double)1 / 3, (double)2 / 3, 1}, {0, 0, (double)1 / 3, 1}, {0, 0, 0, 1} };

        return t;
    }

    static Mat BeizerBasis()
    {
        Mat t = { {1, -3, 3, -1}, {0, 3, -6, 3}, {0, 0, 3, -3}, {0, 0, 0, 1} };
        return t;
    }

    static Mat HermiteBasis()
    {
        //Mat t(4, 4);
        //t[0][0] = 2;
        //t[0][1] = -3;
        //t[0][3] = 1;

        //t[1][0] = -2;
        //t[1][1] = 3;

        //t[2][1] = 1;
        //t[2][2] = -2;
        //t[2][3] = 1;

        //t[3][1] = 1;
        //t[3][2] = -1;

        Mat t = { {1, 0, -3, 2}, {0, 1, -2, 1}, {0, 0, 3, -2}, {0, 0, -1, 1} };

        return t;
    }

};

Mat operator * (Mat a, Mat b)
{
    return a.Mult(b);
}

double _DEGREES2RADIANS(long double a)
{
    return a * DEGREES2RADIANS;
}

double _RADIANS2DEGREES(long double r)
{
    return r * RADIANS2DEGREES;
}

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
    if (isNeg)
    {
        c = a.substr(0, accuracy + 3);
    }
    else
    {
        c = a.substr(0, accuracy + 2);
    }
    return stod(c);
}

double deriv2angle(double d)
{
    return _RADIANS2DEGREES(atan(d));
}

struct Waypoint
{
    double X, Y, Angle;
    Waypoint(double x, double y, double angleInDegrees) : X(x), Y(y), Angle(angleInDegrees) {}
    Waypoint() : Waypoint(0, 0, 0) {}
};

inline bool operator==(const Waypoint& lhs, const Waypoint& rhs)
{
    if (lhs.X == rhs.X && lhs.Y == rhs.Y && lhs.Angle == rhs.Angle)
    {
        return true;
    }
    else
    {
        return false;
    }
}

typedef std::vector<Waypoint> Path;

struct StandardCubicFunction
{
    double A, B, C, D;
    double translateX, translateY;
    double valueAt(const double x) const
    {
        return (A * std::pow(x - translateX, 3)) + (B * std::pow(x - translateX, 2)) + (C * (x - translateX)) + (D + translateY);
    }
    double derivAt(const double x) const
    {
        return (3 * A * std::pow(x - translateX, 2)) + (2 * B * (x - translateX)) + (C);
    }
    StandardCubicFunction(double a, double b, double c, double d, double tx = 0, double ty = 0) : A(a), B(b), C(c), D(d), translateX(tx), translateY(ty) {}
    StandardCubicFunction() : StandardCubicFunction(0, 0, 0, 0) {}
};

struct Spline;
Spline HermiteFinder(Waypoint PointOne, Waypoint PointTwo);
struct Spline
{
    Waypoint point1, point2;
    StandardCubicFunction function;
    //accounts for tx & ty
    Waypoint realPoint1()
    {
        return Waypoint(point1.X - function.translateX, point1.Y + function.translateY, point1.Angle);
    }

    Waypoint realPoint2()
    {
        return Waypoint(point2.X - function.translateX, point2.Y + function.translateY, point2.Angle);
    }
    Spline(Waypoint p1, Waypoint p2, double A, double B, double C, double D, double tx = 0, double ty = 0) : point1(p1), point2(p2), function(StandardCubicFunction(A, B, C, D, tx, ty)) {}
    Spline(Waypoint p1, Waypoint p2, StandardCubicFunction func) : point1(p1), point2(p2), function(func) {}
    Spline(Waypoint p1, Waypoint p2) : Spline(HermiteFinder(p1, p2)) {}
    Spline() : Spline(Waypoint(0, 0, 0), Waypoint(0, 0, 0), 0, 0, 0, 0) {}
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
    if (PointOne == PointTwo)
    {
        return Spline();
    }
    double* A0 = new double(0), * A1 = new double(0), * A2 = new double(0), * A3 = new double(0);
    hermite_cubic_to_power_cubic(PointOne.X, PointOne.Y, Angle2Deriv(PointOne.Angle), PointTwo.X, PointTwo.Y, Angle2Deriv(PointTwo.Angle), A0, A1, A2, A3);
    return Spline(PointOne, PointTwo, truncate(*A3, 10), truncate(*A2, 10), truncate(*A1, 10), truncate(*A0, 10));
}

double distance(Waypoint p1, Waypoint p2)
{
    return std::sqrt(std::pow(p2.X - p1.X, 2) + std::pow(p2.Y - p1.Y, 2));
}

double distance(double x1, double y1, double x2, double y2)
{
    return distance(Waypoint(x1, y1, 0), Waypoint(x2, y2, 0));
}

//uses arc length formula to find distance
double ArcLengthDistance(Spline Function, size_t Accuracy = 500)
{
    double x1 = Function.point1.X;
    double x2 = Function.point2.X;
    ;

    auto f = [Function](double x) {
        return Function.function.valueAt(x);
    };

    double arcLength = 0;
    double h = (x2 - x1) / Accuracy;
    for (size_t i = 0; i <= Accuracy; i++)
    {
        // arcLength += (h/2) * ((f((x1+i)*h)) + (f(x1+(i*h))));
        // std::cout << "At X Value " << x1 + (h * i) << " Y equals: " << f(h * i) << "\n";
        auto nx1 = x1 + (h * i);
        auto ny1 = f(x1 + (h * i));
        auto nx2 = x1 + (h * (i + 1));
        auto ny2 = f(x1 + h * (i + 1));
        auto d = distance(nx1, ny1, nx2, ny2);
        arcLength += d;
    }
    return std::abs(arcLength);
}

double ArcLengthDistance(Spline theSplineFunction, double lowerLimit, double upperLimit, size_t Accuracy = 500)
{
    auto p1 = Waypoint(lowerLimit, theSplineFunction.function.valueAt(lowerLimit), std::atan(theSplineFunction.function.derivAt(lowerLimit)) * RADIANS2DEGREES);
    auto p2 = Waypoint(upperLimit, theSplineFunction.function.valueAt(upperLimit), std::atan(theSplineFunction.function.derivAt(upperLimit)) * RADIANS2DEGREES);
    auto s = HermiteFinder(p1, p2);
    return ArcLengthDistance(
        HermiteFinder(p1, p2), Accuracy);
}

//Parametrizes the cubic function to use arc length and returns the x value arcLength away
double ArcLengthToXValue(Spline spline, double x, double arcLength, const double increment = .01, size_t arcLengthAccuracy = 2500)
{
    bool negative = false;
    if (arcLength < 0)
    {
        negative = true;
    }

    double i = x;
    if (arcLength == 0)
    {
        return 0;
    }

    if (negative)
    {
        while (ArcLengthDistance(spline, i, x, arcLengthAccuracy) <= -arcLength)
        {
            i -= increment;
        }
    }
    else
    {
        while (ArcLengthDistance(spline, x, i, arcLengthAccuracy) <= arcLength)
        {
            i += increment;
        }
    }

    return i;
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
    //TIME STARTS AT 0 THAT'S WHY FUNCTION STARTS at 0, 0
    double A, B, C, D, Time;
    Time = TimeGivenSFJ(Function, Jerk);
    Spline XFunction = HermiteFinder(Waypoint(0, 0, Function.point1.Angle), Waypoint(Time, Function.point2.X - Function.point1.X, Function.point2.Angle));
    Spline YFunction = HermiteFinder(Waypoint(0, 0, Function.point1.Angle), Waypoint(Time, Function.point2.Y - Function.point1.Y, Function.point2.Angle));
    A = Function.function.A;
    B = Function.function.B;
    C = Function.function.C;
    D = Function.function.D;
    return Segment(Function, XFunction, YFunction, A, B, C, D, Time);
}

//these values determine how many segments are created.
//USING_DISTANCE_PARTION determines if a new segment is created based off distance or a equal division
//TANK_CONFIG_PARTION_VALUE represents the distance or how many divisions to create
bool USING_DISTANCE_PARTION = true;
size_t TANK_CONFIG_PARTION_VALUE = .25;

Curve curveGenerator(Path path)
{
    Curve ReturnSpline;
    size_t NumberOfFunctions = path.size() - 1;
    for (int i = 0; i < NumberOfFunctions; i++)
    {
        auto Temp = HermiteFinder(path[i], path[i + 1]);
        ReturnSpline.push_back(Temp);
    }
    return ReturnSpline;
}

typedef std::vector<Segment> Trajectory;

class TankConfig
{
    Curve curve;

    Trajectory leftTrajectory;
    Trajectory rightTrajectory;

    double currentTime = 0.0;

    friend void createDesmosGraph(TankConfig config, std::string fileName, std::string htmlFileLocation);

public:
    enum GlobalType
    {
        VELOCITY,
        ACCELERATION,
        JERK
    };

    void run(double (*RPMFn)())
    {
        std::cout << RPMFn() << "\n";
    }

    template <typename Fn, typename Src>
    void run(Fn&& RPMFn, Src&& src)
    {
        auto f = std::bind(RPMFn, src);
        std::cout << f() << "\n";
    }

    //this assumes PIDCalcFN takes in order of current, target; assumes SetFn is within the src object
    template <typename Fn, typename Src, typename PIDFn, typename PIDSrc, typename SetFN>
    void run(Fn&& RPMFn, Src&& src, PIDFn&& PIDCalcFn, PIDSrc&& PIDController, SetFN&& setFn)
    {
        auto target = 300;
        auto getRPM = std::bind(RPMFn, src);
        auto calc = std::bind(PIDCalcFn, PIDController, std::placeholders::_1, target);
        auto set = std::bind(setFn, src, std::placeholders::_1);
        auto s = getRPM();
        auto e = calc(s);
        std::cout << "Speed is: " << s << "\n";
        std::cout << "New voltage is: " << e << "\n";
        set(e);
    }

    void testTrajectory()
    {
        std::cout << "working" << "\n";
        for (size_t i = 0; i < leftTrajectory.size(); i++)
        {
            std::cout << "Left Segment " << i << "'s Values: " << leftTrajectory[i].A << " " << leftTrajectory[i].B << " " << leftTrajectory[i].C << " " << leftTrajectory[i].D << " " << leftTrajectory[i].spline.point1.X << " " << leftTrajectory[i].spline.point1.Y << "\n";
        }
        for (size_t i = 0; i < rightTrajectory.size(); i++)
        {
            std::cout << "Right Segment " << i << "'s Values: " << rightTrajectory[i].A << " " << rightTrajectory[i].B << " " << rightTrajectory[i].C << " " << rightTrajectory[i].D << "\n";
        }
    }

    static const int parts = 8;
    //this is just wrong, need to just get points away from function and make that left and right
    TankConfig(Curve curve, double widthBetweenWheels, double jerk)
    {
        TankConfig::curve = curve;
        leftTrajectory = Trajectory();
        rightTrajectory = Trajectory();
        double radius = widthBetweenWheels / 2;

        Path leftPath, rightPath;

        for (size_t i = 0; i < curve.size(); i++)
        {
            auto s = curve[i];
            auto p1 = s.point1, p2 = s.point2;
            Mat ctrl(4, 1);
            ctrl[0][0] = p1.X;
            ctrl[0][1] = Angle2Deriv(p1.Angle);
            ctrl[0][2] = p2.X;
            ctrl[0][3] = Angle2Deriv(p2.Angle);

            Mat conversion = Mat::InverseBeizerBasis() * Mat::HermiteBasis();
            Mat bezierCtrl = conversion * ctrl;

            Mat t(1, 4, 1);
            t.setPowers({ {1}, {1}, {2}, {3} });
            t.setVariables({ {0}, {1}, {1}, {1} });

            //actual number below
            auto result = t.input(.5) * Mat::BeizerBasis() * bezierCtrl;

        }

        auto leftCurve = curveGenerator(leftPath);
        auto rightCurve = curveGenerator(rightPath);

        for (size_t i = 0; i < leftCurve.size(); i++)
        {
            leftTrajectory.push_back(GenerateSegment(leftCurve[i], jerk));
            rightTrajectory.push_back(GenerateSegment(rightCurve[i], jerk));
        }
    }

    double globalAt(double time, GlobalType type, bool left)
    {
        double t = 0;
        size_t i = 0;
        if (left)
        {
            while (t > time)
            {
                time += leftTrajectory[i].time;
                i++;
            }
            switch (type)
            {
            case GlobalType::VELOCITY:
                return leftTrajectory[i].Velocity(t - time);
            case GlobalType::ACCELERATION:
                return leftTrajectory[i].Acceleration(t - time);
            case GlobalType::JERK:
                return leftTrajectory[i].Jerk(t - time);
            }
        }
        else
        {
            while (t > time)
            {
                time += rightTrajectory[i].time;
                i++;
            }
            switch (type)
            {
            case GlobalType::VELOCITY:
                return rightTrajectory[i].Velocity(t - time);
            case GlobalType::ACCELERATION:
                return rightTrajectory[i].Acceleration(t - time);
            case GlobalType::JERK:
                return rightTrajectory[i].Jerk(t - time);
            }
        }
    }
};

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
    file << "<!DOCTYPE html><html><head></head><style>html, body{margin: 0px; width: 100%; height: 100%;}</style><body>";
    file << u8R"(<script src="https://www.desmos.com/api/v1.5/calculator.js?apiKey=dcb31709b452b1cf9dc26972add0fda6"></script>)"s;

    file << u8R"(<div id="calculator" style="width: calc(100% - 2px); height: calc(100% - 2px);"></div>)"s;

    file << "<script>\n";
    file << u8R"(var elt = document.getElementById('calculator');)"s << "\n"
        << "var calculator = Desmos.GraphingCalculator(elt);";

    size_t precision = 10;

    for (size_t i = 0; i < config.leftTrajectory.size(); i++)
    {
        auto s = config.leftTrajectory[i];
        file << u8R"(calculator.setExpression({id: 'leftGraph)"s;
        file << i;
        file << u8R"(', color: Desmos.Colors.RED, latex: 'y=)"s;
        file << fixed << setprecision(precision) << s.A << "(x - " << s.spline.function.translateX << ")^3 + " << setprecision(precision) << fixed << s.B << "(x -" << s.spline.function.translateX << ")^2 + " << setprecision(precision) << fixed << s.C << "(x -" << s.spline.function.translateX << ") + " << setprecision(precision) << fixed << (s.D + s.spline.function.translateY);
        file << "  \\\\left\\\\{" << setprecision(6) << min(s.spline.realPoint1().X, s.spline.realPoint2().X) << "\\\\le x \\\\le" << setprecision(6) << max(s.spline.realPoint2().X, s.spline.realPoint2().X) << "\\\\right\\\\}";
        file << u8R"('});)"s;
    }

    for (size_t i = 0; i < config.rightTrajectory.size(); i++)
    {
        auto s = config.rightTrajectory[i];
        file << u8R"(calculator.setExpression({id: 'rightGraph)"s;
        file << i;
        file << u8R"(', color: Desmos.Colors.BLUE, latex: 'y=)"s;
        file << fixed << setprecision(precision) << s.A << "(x - " << s.spline.function.translateX << ")^3 + " << setprecision(precision) << fixed << s.B << "(x -" << s.spline.function.translateX << ")^2 + " << setprecision(precision) << fixed << s.C << "(x -" << s.spline.function.translateX << ") + " << setprecision(precision) << fixed << (s.D + s.spline.function.translateY);
        file << "  \\\\left\\\\{" << setprecision(6) << min(s.spline.realPoint1().X, s.spline.realPoint2().X) << "\\\\le x \\\\le" << setprecision(6) << max(s.spline.realPoint2().X, s.spline.realPoint2().X) << "\\\\right\\\\}";
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

    file << u8R"(<div id="calculator" style="width: calc(100% - 2px); height: calc(100% - 2px);"></div>)"s;

    file << "<script>\n";
    file << u8R"(var elt = document.getElementById('calculator');)"s << "\n"
        << "var calculator = Desmos.GraphingCalculator(elt);";

    size_t precision = 10;

    for (size_t i = 0; i < curve.size(); i++)
    {
        file << u8R"(calculator.setExpression({id: 'test)"s << i << "'";
        file << u8R"(, color: Desmos.Colors.BLACK, latex: 'y=)"s;

        file << setprecision(precision) << fixed << curve[i].function.A << "x^3 + " << setprecision(precision) << fixed
            << curve[i].function.B << "x^2 + " << setprecision(precision) << fixed
            << curve[i].function.C << "x + " << setprecision(precision) << fixed << curve[i].function.D;

        file << "  \\\\left\\\\{" << setprecision(6) << min(curve[i].point1.X, curve[i].point2.X)
            << "\\\\le  x \\\\le" << setprecision(6) << max(curve[i].point1.X, curve[i].point2.X) << "\\\\right\\\\} ";

        file << u8R"('});)"s;
    }
    file << "</script>";

    file.close();
}

template <typename SpeedController>
class TankDrive
{
    TankConfig* config;
    double time;
    SpeedController& masterLeft, & masterRight;
    std::function<double(double, double)> calc;
    std::function<double()> getLeft;
    std::function<double()> getRight;
    std::function<void(double)> setLeft;
    std::function<void(double)> setRight;

public:
    TankDrive(SpeedController& leftMaster, SpeedController& rightMaster, TankConfig* con) : masterLeft(leftMaster), masterRight(rightMaster), config(con), time(0)
    {
    }

    void setTime(double t)
    {
        time = t;
    }

    //assumes actual value, target
    template <typename PIDFn, typename Src>
    void setUpPID(PIDFn&& fn, Src&& src)
    {
        calc = std::bind(fn, src, std::placeholders::_1, std::placeholders::_2);
    }

    template <typename getFn, typename Src>
    void setUpGet(getFn&& fn, Src&& leftSrc, Src&& rightSrc)
    {
        getLeft = std::bind(fn, leftSrc);
        getRight = std::bind(fn, rightSrc);
    }

    template <typename setFn, typename Src>
    void setUpSet(setFn&& fn, Src&& leftSrc, Src&& rightSrc)
    {
        setLeft = std::bind(fn, leftSrc, std::placeholders::_1);
        setRight = std::bind(fn, rightSrc, std::placeholders::_1);
    }

    void run()
    {
        auto i1 = calc(getLeft(), config->globalAt(time, TankConfig::GlobalType::VELOCITY, true));
        auto i2 = calc(getRight(), config->globalAt(time, TankConfig::GlobalType::VELOCITY, true));
        std::cout << "Set Voltage are: " << i1 << " " << i2 << "\n";
        setLeft(i1);
        setRight(i2);
        time += 0.02;
    }
};