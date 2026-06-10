#ifndef RIEMANN_H_
#define RIEMANN_H_

#include <vector>

void exact_riemann (std::vector<double>& x, std::vector<double>& p,
                    std::vector<double>& rho, std::vector<double>& v, double t,
                    double pl, double rhol, double vl,
                    double pr, double rhor, double vr);

#endif
