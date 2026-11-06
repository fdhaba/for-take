// ============================================================================
// Assignment 2 - Numerical Methods: Root-Finding Algorithms
// Course : Computational Methods
// Methods: Bisection, False Position, Newton-Raphson, Secant, Fixed-Point
// ============================================================================
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <string>
#include <limits>
 
using namespace std;
 

const string STUDENT_NAME = "Tekalign Birhanu";
const string STUDENT_ID   = "RU2091/17";
// --------------------------------------------------
 
// ---------------------------- Test functions -------------------------------
// Each entry: f(x), f'(x), g(x) for fixed point (x = g(x)), default inputs.
struct Problem {
    string name;
    double (*f)(double);
    double (*df)(double);
    double (*g)(double);
    double a, b, x0;   // default bracket and initial guess
};
 
double f1(double x)  { return x*x*x - 2*x - 5; }
double df1(double x) { return 3*x*x - 2; }
double g1(double x)  { return cbrt(2*x + 5); }
 
double f2(double x)  { return cos(x) - x; }
double df2(double x) { return -sin(x) - 1; }
double g2(double x)  { return cos(x); }
 
double f3(double x)  { return exp(-x) - x; }
double df3(double x) { return -exp(-x) - 1; }
double g3(double x)  { return exp(-x); }
 
vector<Problem> problems = {
    {"f(x) = x^3 - 2x - 5   [g(x) = cbrt(2x+5)]", f1, df1, g1, 2.0, 3.0, 2.0},
    {"f(x) = cos(x) - x     [g(x) = cos(x)]",     f2, df2, g2, 0.0, 1.0, 0.5},
    {"f(x) = e^(-x) - x     [g(x) = e^(-x)]",     f3, df3, g3, 0.0, 1.0, 0.5}
};
 
// ------------------------------ Result type --------------------------------
struct Row { int iter; double xr, fxr, ea; bool hasEa; };
 
struct Result {
    string method;
    vector<Row> rows;
    bool converged = false;
    string note;
};
 
double relErrPct(double xnew, double xold) {
    if (xnew == 0.0) return fabs(xnew - xold) * 100.0;
    return fabs((xnew - xold) / xnew) * 100.0;
}
 
// ------------------------------- Methods -----------------------------------
Result bisection(const Problem& p, double a, double b, double tol, int maxIt) {
    Result r; r.method = "Bisection";
    if (p.f(a) * p.f(b) >= 0) { r.note = "Invalid bracket: f(a)*f(b) must be < 0."; return r; }
    double xr = a, xold = a;
    for (int i = 1; i <= maxIt; ++i) {
        xold = xr;
        xr = (a + b) / 2.0;
        double fx = p.f(xr);
        bool has = (i > 1);
        double ea = has ? relErrPct(xr, xold) : 0.0;
        r.rows.push_back({i, xr, fx, ea, has});
        if (fx == 0.0 || (has && ea < tol)) { r.converged = true; break; }
        if (p.f(a) * fx < 0) b = xr; else a = xr;
    }
    if (!r.converged) r.note = "Max iterations reached without meeting tolerance.";
    return r;
}
 
Result falsePosition(const Problem& p, double a, double b, double tol, int maxIt) {
    Result r; r.method = "False Position";
    if (p.f(a) * p.f(b) >= 0) { r.note = "Invalid bracket: f(a)*f(b) must be < 0."; return r; }
    double xr = a, xold = a;
    for (int i = 1; i <= maxIt; ++i) {
        xold = xr;
        double fa = p.f(a), fb = p.f(b);
        xr = b - fb * (a - b) / (fa - fb);
        double fx = p.f(xr);
        bool has = (i > 1);
        double ea = has ? relErrPct(xr, xold) : 0.0;
        r.rows.push_back({i, xr, fx, ea, has});
        if (fx == 0.0 || (has && ea < tol)) { r.converged = true; break; }
        if (fa * fx < 0) b = xr; else a = xr;
    }
    if (!r.converged) r.note = "Max iterations reached without meeting tolerance.";
    return r;
}
 
Result newtonRaphson(const Problem& p, double x0, double tol, int maxIt) {
    Result r; r.method = "Newton-Raphson";
    double x = x0;
    for (int i = 1; i <= maxIt; ++i) {
        double d = p.df(x);
        if (fabs(d) < 1e-14) { r.note = "Derivative is (near) zero; method fails."; return r; }
        double xn = x - p.f(x) / d;
        double ea = relErrPct(xn, x);
        r.rows.push_back({i, xn, p.f(xn), ea, true});
        x = xn;
        if (ea < tol || p.f(xn) == 0.0) { r.converged = true; break; }
    }
    if (!r.converged && r.note.empty()) r.note = "Max iterations reached without meeting tolerance.";
    return r;
}
 
Result secant(const Problem& p, double x0, double x1, double tol, int maxIt) {
    Result r; r.method = "Secant";
    double xm = x0, x = x1;
    for (int i = 1; i <= maxIt; ++i) {
        double fm = p.f(xm), fx = p.f(x);
        if (fabs(fx - fm) < 1e-14) { r.note = "Division by (near) zero in secant formula."; return r; }
        double xn = x - fx * (xm - x) / (fm - fx);
        double ea = relErrPct(xn, x);
        r.rows.push_back({i, xn, p.f(xn), ea, true});
        xm = x; x = xn;
        if (ea < tol || p.f(xn) == 0.0) { r.converged = true; break; }
    }
    if (!r.converged && r.note.empty()) r.note = "Max iterations reached without meeting tolerance.";
    return r;
}
 
Result fixedPoint(const Problem& p, double x0, double tol, int maxIt) {
    Result r; r.method = "Fixed-Point Iteration";
    double x = x0;
    for (int i = 1; i <= maxIt; ++i) {
        double xn = p.g(x);
        if (!isfinite(xn)) { r.note = "Iteration diverged (non-finite value)."; return r; }
        double ea = relErrPct(xn, x);
        r.rows.push_back({i, xn, p.f(xn), ea, true});
        x = xn;
        if (ea < tol) { r.converged = true; break; }
    }
    if (!r.converged && r.note.empty()) r.note = "Max iterations reached without meeting tolerance.";
    return r;
}
 
// ------------------------------- Output ------------------------------------
void printTable(const Result& r) {
    cout << "\n--- " << r.method << " ---\n";
    cout << left << setw(6) << "Iter" << setw(18) << "x_r" << setw(18) << "f(x_r)"
         << setw(16) << "% Ea" << "\n";
    cout << string(58, '-') << "\n";
    for (const auto& row : r.rows) {
        cout << left << setw(6) << row.iter
             << setw(18) << fixed << setprecision(8) << row.xr
             << setw(18) << scientific << setprecision(4) << row.fxr;
        if (row.hasEa) cout << setw(16) << scientific << setprecision(4) << row.ea;
        else           cout << setw(16) << "---";
        cout << "\n";
    }
    cout << fixed;
    if (!r.rows.empty())
        cout << "Result: root ~ " << setprecision(8) << r.rows.back().xr
             << " after " << r.rows.size() << " iteration(s)";
    cout << (r.converged ? "  [CONVERGED]" : "  [DID NOT CONVERGE]") << "\n";
    if (!r.note.empty()) cout << "Note: " << r.note << "\n";
}
 
void printSummary(const vector<Result>& all, double tol) {
    cout << "\n====================== COMPARISON SUMMARY ======================\n";
    cout << "Stopping tolerance: % Ea < " << tol << "\n";
    cout << left << setw(24) << "Method" << setw(12) << "Iterations" << setw(18) << "Root estimate"
         << setw(14) << "Status" << "\n";
    cout << string(68, '-') << "\n";
    for (const auto& r : all) {
        cout << left << setw(24) << r.method << setw(12) << r.rows.size();
        if (!r.rows.empty()) cout << setw(18) << fixed << setprecision(8) << r.rows.back().xr;
        else cout << setw(18) << "n/a";
        cout << setw(14) << (r.converged ? "Converged" : "Failed") << "\n";
    }
}
 
template <typename T>
T ask(const string& prompt, T def) {
    cout << prompt << " [default " << def << "]: ";
    string line; getline(cin, line);
    if (line.empty()) return def;
    try {
        if constexpr (is_same<T, int>::value) return stoi(line);
        else return stod(line);
    } catch (...) { return def; }
}
 
int main() {
    cout << "==========================================================\n";
    cout << "  ROOT-FINDING ALGORITHMS - Computational Methods, Assignment 2\n";
    cout << "  Student Name : " << STUDENT_NAME << "\n";
    cout << "  Student ID   : " << STUDENT_ID << "\n";
    cout << "==========================================================\n\n";
 
    cout << "Select the equation f(x) = 0:\n";
    for (size_t i = 0; i < problems.size(); ++i)
        cout << "  " << i + 1 << ") " << problems[i].name << "\n";
    int choice = ask<int>("Choice (1-3)", 1);
    if (choice < 1 || choice > (int)problems.size()) choice = 1;
    const Problem& p = problems[choice - 1];
 
    cout << "\nEnter inputs (press Enter to accept the default shown):\n";
    double a   = ask<double>("Bracket lower bound a", p.a);
    double b   = ask<double>("Bracket upper bound b", p.b);
    double x0  = ask<double>("Initial guess x0 (Newton, Fixed-Point)", p.x0);
    double tol = ask<double>("Tolerance, % relative error", 0.01);
    int maxIt  = ask<int>("Maximum iterations", 100);
    cout << "(Secant uses a and b as its two starting points.)\n";
 
    vector<Result> all;
    all.push_back(bisection(p, a, b, tol, maxIt));
    all.push_back(falsePosition(p, a, b, tol, maxIt));
    all.push_back(newtonRaphson(p, x0, tol, maxIt));
    all.push_back(secant(p, a, b, tol, maxIt));
    all.push_back(fixedPoint(p, x0, tol, maxIt));
 
    cout << "\nEquation: " << p.name << "\n";
    cout << "Tolerance: " << tol << " %   Max iterations: " << maxIt << "\n";
    for (const auto& r : all) printTable(r);
    printSummary(all, tol);
 
    cout << "\nStudent: " << STUDENT_NAME << "  |  ID: " << STUDENT_ID << "\n";
    cout << "\nPress Enter to exit...";
    string dummy; getline(cin, dummy);
    return 0;
}