#include "AutoPilot.hpp"
#include <fstream>
#include <iomanip>

using namespace std::string_literals;
using namespace std;

int main() {
    cout << "Tesing AP\n";
    cout.precision(6);
    Path path = {Waypoint(0, 0, 30), Waypoint(1, 1, 45)};
    auto curve = curveGenerator(path);
    TankConfig drive(curve, 2, 2);
    drive.testTrajectory();
    ofstream file;
    file.open("C:\\Users\\Lance\\Documents\\AutoPilot-Library\\data.html");
    file << "<!DOCTYPE html><html><head></head><body>";
    file << u8R"(<script src="https://www.desmos.com/api/v1.5/calculator.js?apiKey=dcb31709b452b1cf9dc26972add0fda6"></script>)"s;


    file << u8R"(<div id="calculator" style="width: 1200px; height: 800px;"></div>)"s;
    
    file << "<script>\n";
    file << u8R"(var elt = document.getElementById('calculator');)"s << "\n" << "var calculator = Desmos.GraphingCalculator(elt);";

    for (size_t i = 0; i < drive.leftTrajectory.size(); i++)
    {
        auto s = drive.leftTrajectory[i];
        file << u8R"(calculator.setExpression({id: 'leftGraph)"s;
        file << i;
        file << u8R"(', latex: 'y=)"s;
        file << fixed << setprecision(20) << s.A << "x^3 + " << setprecision(20) << fixed << s.B << "x^2 + " << setprecision(20) << fixed << s.C << "x + " << setprecision(20) << fixed << s.D;
        file << " \\\\left\\\\{" << setprecision(6) << min(s.spline.point1.X, s.spline.point2.X) << "<=x<=" << setprecision(6) << max(s.spline.point1.X, s.spline.point2.X) << ":x\\\\right\\\\}";
        file << u8R"('});)"s;
    }

    for (size_t i = 0; i < drive.rightTrajectory.size(); i++)
    {
        auto s = drive.rightTrajectory[i];
        file << u8R"(calculator.setExpression({id: 'rightGraph)"s;
        file << i;
        file << u8R"(', latex: 'y=)"s;
        file << setprecision(20) << fixed << s.A << "x^3 + " << setprecision(20) << fixed << s.B << "x^2 + " << setprecision(20) << fixed << s.C << "x + " << setprecision(20) << fixed << s.D;
        file << " \\\\left\\\\{" << setprecision(6) << min(s.spline.point1.X, s.spline.point2.X) << "<=x<=" << setprecision(6) << max(s.spline.point1.X, s.spline.point2.X) << ":x\\\\right\\\\}";
        file << u8R"('});)"s;
    }

    auto spline = HermiteFinder(path[0], path[1]);
    file << u8R"(calculator.setExpression({id: 'test')"s;
    file << u8R"(, latex: 'y=)"s;
    file << setprecision(20) << fixed << spline.function.A << "x^3 + " << setprecision(20) << fixed << spline.function.B << "x^2 + " << setprecision(20) << fixed << spline.function.C << "x + " << setprecision(20) << fixed << spline.function.D;
    file << " \\\\left\\\\{" << setprecision(6) << min(spline.point1.X, spline.point2.X) << "<=x<=" << setprecision(6) << max(spline.point1.X, spline.point2.X) << ":x\\\\right\\\\}";
    file << u8R"('});)"s;

    file << "</script>";
    // file << u8R"(
    // <script>
    //     var elt = document.getElementById('calculator');
    //     var calculator = Desmos.GraphingCalculator(elt);
    //     calculator.setExpression({ id: 'graph1', latex: 'y=x^2' });
    // </script>
    // )"s;

    file.close();
    return 0;
} 