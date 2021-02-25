#include "AutoPilot.hpp"

using namespace std;

const string desktop = "C:\\Users\\Lance\\Documents\\AutoPilot-Library\\";
const string laptop = "C:\\Users\\Lance\\Documents\\AutoPilot-Library\\";

int main(int argc, char *argv[]) {


        string PATH = getenv("PROGRAMPATH");
        PATH.append("\\");
        cout << PATH << "\n";

        cout << "Testing AP\n";
        cout.precision(10);
        Path path = {Waypoint(0, 0, 0), Waypoint(1, 2, 0), Waypoint(3, 4, 0), Waypoint(5, 0, 0)};
        auto curve = curveGenerator(path);
        TankConfig drive(curve, 2, 1);
        // drive.testTrajectory();
        createDesmosGraph(drive, "", PATH + "graph.html");
        // Spline XSpline = HermiteFinder(Waypoint( 0, path[0].X, 45), Waypoint( path[1].X - path[0].X, path[1].X, 45));
        // Spline YSpline = HermiteFinder(Waypoint( 0, path[0].Y, path[0].Angle ), Waypoint( path[1].X - path[0].X, path[1].Y, path[1].Angle ));
        // createDesmosGraph(Curve{XSpline}, "", desktop + "graph1.html");
        // createDesmosGraph(Curve{YSpline}, "", desktop + "graph2.html");
        // cout << "ArcLengthTest: " << ArcLengthDistance(curve[0], 1, -1) << "\n";
        // cout << "ArcLengthToXValue Test: " << ArcLengthToXValue(curve[0], 1, -.25) << "\n";
        // cout << "Distance test: " << distance(0, 0, 2, 3    ) << "\n";
        return 0;
}

//maybe do a bit of automation for increased smoothness. get linear slope between two slopes and like use that instead of user input if not given
