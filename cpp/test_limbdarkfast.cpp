// Check the C++ version against the shared test vectors.
//
//     ./test_limbdarkfast [path/to/tests/vectors]
//
// Exits non-zero if any point differs from the reference by more than TOL.

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "limbdarkfast.hpp"

static const double TOL = 1e-12;

static std::vector<std::vector<double> > read_rows(const std::string& path) {
    std::ifstream in(path.c_str());
    if (!in) throw std::runtime_error("cannot open " + path);
    std::vector<std::vector<double> > rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::vector<double> row;
        double v;
        while (ss >> v) row.push_back(v);
        rows.push_back(row);
    }
    return rows;
}

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? argv[1] : "../tests/vectors";
    const char* laws[3] = {"hybrid2", "hybrid4", "hybrid5"};
    const char* kinds[2] = {"sphere", "oblate"};
    int failures = 0;
    for (int ik = 0; ik < 2; ++ik) {
        for (int il = 0; il < 3; ++il) {
            const std::string kind = kinds[ik];
            const limbdarkfast::Law law = limbdarkfast::law_from_name(laws[il]);
            const int n_w = limbdarkfast::n_weights(law);
            const int n_geo = kind == "sphere" ? 2 : 5;
            const std::vector<std::vector<double> > rows =
                read_rows(dir + "/" + kind + "_" + laws[il] + ".txt");
            double worst = 0.0;
            for (std::size_t i = 0; i < rows.size(); ++i) {
                const std::vector<double>& d = rows[i];
                const double* w = &d[n_geo];
                const double F = kind == "sphere"
                    ? limbdarkfast::hybrid_flux(d[0], d[1], w, law)
                    : limbdarkfast::hybrid_flux_oblate(d[0], d[1], d[2], d[3], d[4], w, law);
                const double err = std::fabs(F - d[n_geo + n_w]);
                if (!(err <= worst)) worst = err;      // NaN counts as worst
            }
            const bool ok = worst <= TOL;
            failures += !ok;
            std::printf("%-6s %s: %5zu points, max |dF| = %.2e  %s\n", kind.c_str(),
                        laws[il], rows.size(), worst, ok ? "ok" : "FAIL");
        }
    }
    // every shape is 1 at the limb and 0 at the centre
    for (int il = 0; il < 3; ++il) {
        const limbdarkfast::Law law = limbdarkfast::law_from_name(laws[il]);
        for (int j = 0; j < limbdarkfast::n_weights(law); ++j) {
            double w[5] = {0, 0, 0, 0, 0};
            w[j] = 1.0;
            if (std::fabs(limbdarkfast::intensity(0.0, w, law)) > 1e-12 ||
                std::fabs(limbdarkfast::intensity(1.0, w, law) - 1.0) > 1e-12) {
                std::printf("shape %d of %s is not 1 at the limb and 0 at the centre\n",
                            j + 1, laws[il]);
                ++failures;
            }
        }
    }
    std::printf("law shapes: %s\n", failures ? "see above" : "ok");
    return failures ? 1 : 0;
}
