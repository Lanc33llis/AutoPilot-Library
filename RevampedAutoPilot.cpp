#include "AutoPilot.hpp"

using namespace std;

const string desktop = "C:\\Users\\Lance\\Documents\\AutoPilot-Library\\";
const string laptop = "C:\\Users\\Lance\\Documents\\AutoPilot-Library\\";


struct testSC
{
    double rpm;
    double voltage;
    testSC() : rpm(.0), voltage(.0) {};
    double getVelocity()
    {
        rpm += (voltage*2);
        return rpm*4;
    }
    void set(double i)
    {
        voltage = i;
    }
};

struct testPIDC
{
    double p, i, d;
    testPIDC()
    {
        p = .15, i = 0, d = 0;
    }
    double calculate(double current, double target)
    {
        return p * (target - current);
    }
};

double testRPM()
{
    static double i = 0;
    i++;
    return i*4;
}


int main() {
    cout << "Tesing AP\n";
    cout.precision(10);
    Path path = {Waypoint(0, 0, 0), Waypoint(1, 2, 0), Waypoint(3, 4, 0), Waypoint(5, 0, 0)};
    auto curve = curveGenerator(path);
    TankConfig drive(curve, 2, 1);
    testSC test;
    testPIDC pid;
    testSC m1, m2;
    TankDrive<testSC> t1(m1, m2, &drive);
    t1.setUpPID(&testPIDC::calculate, pid);
    t1.setUpGet(&testSC::getVelocity, m1, m2);
    t1.setUpSet(&testSC::set, m1, m2);
    for (size_t i = 0; i < 50; i++)
    {
        // std::cout << "Cycle: " << i + 1 << "\n";
        // drive.run(&testSC::getVelocity, &test, &testPIDC::calculate, pid, &testSC::set);
        t1.run();
    }

    // drive.testTrajectory();
    createDesmosGraph(drive, "", desktop + "graph.html");
    // Spline XSpline = HermiteFinder(Waypoint( 0, path[0].X, 45), Waypoint( path[1].X - path[0].X, path[1].X, 45));
    // Spline YSpline = HermiteFinder(Waypoint( , path[0].Y, path[0].Angle ), Waypoint( path[1].X - path[0].X, path[1].Y, path[1].Angle ));
    // createDesmosGraph(Curve{XSpline}, "", desktop + "graph1.html");
    // createDesmosGraph(Curve{YSpline}, "", desktop + "graph2.html");
    // cout << "ArcLengthTest: " << ArcLengthDistance(curve[0], 1, -1) << "\n";
    // cout << "ArcLengthToXValue Test: " << ArcLengthToXValue(curve[0], 1, -.25) << "\n";
    // cout << "Distance test: " << distance(0, 0, 2, 3    ) << "\n";
    return 0;
}

//maybe do a bit of automation for increased smoothness. get linear slope between two slopes and like use that instead of user input if not given
