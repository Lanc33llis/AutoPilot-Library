#include "AutoPilot.hpp"

using namespace std;

const string desktop = "D:\\Documents\\VSAutoPilot\\AutoPilot\\AutoPilot\\";
const string laptop = "C:\\Users\\Lance\\Documents\\AutoPilot-Library\\";


struct testSC
{
    double rpm;
    double voltage;
    testSC() : rpm(.0), voltage(.0) {};
    double getVelocity()
    {
        rpm += (voltage * 2);
        return rpm * 4;
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
    return i * 4;
}


int main() {
    //cout << "Tesing AP\n";
    //cout.precision(10);
    Path path = { Waypoint(-1, 0, 0), Waypoint(1, 2, 0), Waypoint(3, 4, 0), Waypoint(5, 0, 0) };
    //auto curve = curveGenerator(path);
    //TankConfig drive(curve, 2, 1);
    //testSC test;
    //testPIDC pid;
    //testSC m1, m2;
    //TankDrive<testSC> t1(m1, m2, &drive);
    //t1.setUpPID(&testPIDC::calculate, pid);
    //t1.setUpGet(&testSC::getVelocity, m1, m2);
    //t1.setUpSet(&testSC::set, m1, m2);
    //for (size_t i = 0; i < 50; i++)
    //{
    //    // std::cout << "Cycle: " << i + 1 << "\n";
    //    // drive.run(&testSC::getVelocity, &test, &testPIDC::calculate, pid, &testSC::set);
    //    t1.run();
    //}

    //// drive.testTrajectory();
    //createDesmosGraph(drive, "", desktop + "graph.html");

    //return 0;

    auto p1 = path[0], p2 = path[1];
    auto
        x0 = p1.X, x1 = p2.X,
        y0 = p1.Y, y1 = p2.Y,
        d0 = Angle2Deriv(p1.Angle), d1 = Angle2Deriv(p2.Angle);
    auto nx1 = x0, nx2 = x0 + (double)1 / 3 * (x1 - x0), nx3 = x1 - (double)1 / 3 * (x1 - x0), nx4 = x1;
    auto ny1 = y0, ny2 = y0 + (double)1 / 3 * d0 * (y1 - y0), ny3 = y1 + (double)1 / 3 * d1 * (y1 - y0), ny4 = y1;

    Mat yctrl = { {ny1, ny2, ny3, ny4} };
    Mat xctrl = { {nx1, nx2, nx3, nx4} };

    Mat t(1, 4, 1);
    t.setPowers({ {1}, {1}, {2}, {3} });
    t.setVariables({ {0}, {1}, {1}, {1} });
    auto i = t.input(.5);

    auto x = i * Mat::BeizerBasis() * xctrl, y = i * Mat::BeizerBasis() * yctrl;
}

//maybe do a bit of automation for increased smoothness. get linear slope between two slopes and like use that instead of user input if not given
