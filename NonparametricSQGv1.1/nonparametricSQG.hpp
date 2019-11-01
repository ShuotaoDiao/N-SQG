//
//  nonparametricSQG.hpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/23/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#ifndef nonparametricSQG_hpp
#define nonparametricSQG_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <ilcplex/ilocplex.h>
#include <ctime>
#include <stdlib.h>
#include <cassert>
#include <unordered_map>
#include "NonparametricSQG_dataType.hpp"

// test
void CplexConnectionTest(); // test on the cplex connection
void CplexExtremeRayTest(); // test on getting extreme ray
void CplexExtremeRayTest02(); // test on getting extreme ray equality constraint
// obtain dual multiplers of the second stage, given x (first stage decision variable), d, De, Ce, be, Di, Ci, bi
dualMultipliers twoStageLP_secondStageDual(const std::vector<double>& x, const std::vector<double>& d, const std::vector<std::vector<double>>& De, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be, const std::vector<std::vector<double>>& Di, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi);
dualMultipliers twoStageLP_secondStageDual_equality(const std::vector<double>& x, const std::vector<double>& d,const std::vector<std::vector<double>>& Di, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi); // no equality constraint
dualMultipliers twoStageLP_secondStageDual_inequality(const std::vector<double>& x, const std::vector<double>& d,const std::vector<std::vector<double>>& De, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be); // no inequality constraint
// Nonparametric SQG solver for two stage linear programming
std::vector<double> twoStageLP_solver(twoStageParameters parameters,const std::vector<std::vector<std::vector<double>>>& be_database, const std::vector<std::vector<std::vector<double>>>& bi_database, const std::vector<std::vector<double>>& weight_database, double initial_stepsize, int max_iterates, std::vector<double> x_old);
// Nonparametric SQG solver which has output ability
std::vector<double> twoStageLP_solver(twoStageParameters parameters,const std::vector<std::vector<std::vector<double>>>& be_database, const std::vector<std::vector<std::vector<double>>>& bi_database, const std::vector<std::vector<double>>& weight_database, double initial_stepsize, int max_iterates, std::vector<double> x_old, const std::string& resultsOutput_path);
// Nonparametric SQG solver for one stage piecewise linear programming which has output ability
std::vector<double> oneStagePiecewiseLP_solver(oneStageParameters parameters, const oneStagePiecewiseLinearObjectiveDB& parametersDB_objective, double initial_stepsize, int max_iterates, std::vector<double> x_old, const std::string& resultsOutput_path);
// Robust Nonparametric SQG solver which has output ability
std::vector<double> robustTwoStageLP_solver(twoStageParameters parameters,const std::vector<std::vector<std::vector<double>>>& be_database, const std::vector<std::vector<std::vector<double>>>& bi_database, const std::vector<std::vector<double>>& weight_database, int max_iterates, std::vector<double> x_old, double Dx, int m, double M, const std::string& resultsOutput_path);
// Robust Nonparametric SQG solver for one stage piecewise linear programming which has output ability
std::vector<double> robustOneStagePiecewiseLP_solver(oneStageParameters parameters, const oneStagePiecewiseLinearObjectiveDB& parametersDB_objective, int max_iterates, std::vector<double> x_old, double Dx, int m, double M, const std::string& resultsOutput_path);
// for debugging
std::vector<double> robustOneStagePiecewiseLP_solverDEBUG(oneStageParameters parameters, const oneStagePiecewiseLinearObjectiveDB& parametersDB_objective, int max_iterates, std::vector<double> x_old, double Dx, int m, double M, const std::string& resultsOutput_path);
// projection for the first stage in the two stage linear programming
std::vector<double> twoStageLP_projection(const std::vector<double> x, const std::vector<std::vector<double>>& A, const std::vector<double> b, const long& A_rowsize, const long& A_colsize);
// projecting x onto a convex feasible region
std::vector<double> LP_projection(const std::vector<double>& x, const oneStageParameters& parameters, long A_rowsize, long x_size);
/* functions for feasibility cut generation
 */
dualMultipliers twoStageLP_secondStageExtremRay(const std::vector<double>& x, const std::vector<double>& d, const std::vector<std::vector<double>>& De, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be, const std::vector<std::vector<double>>& Di, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi);
dualMultipliers twoStageLP_secondStageExtremRay_inequality(const std::vector<double>& x, const std::vector<double>& d, const std::vector<std::vector<double>>& Di, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi);
dualMultipliers twoStageLP_secondStageExtremRay_equality(const std::vector<double>& x, const std::vector<double>& d, const std::vector<std::vector<double>>& De, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be);
feasibilityCut twoStageLP_feasibilityCutGeneration(const dualMultipliers& extremeRay, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi);
feasibilityCut twoStageLP_feasibilityCutGeneration_inequality(const dualMultipliers& extremeRay, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi);
feasibilityCut twoStageLP_feasibilityCutGeneration_equality(const dualMultipliers& extremeRay, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be);
// Dx finer estimate
double DxEstimate_finer(const std::vector<double>& x, const oneStageParameters& parameters, double lowerBound, double upperBound);
// operator overloading
double operator*(const std::vector<double>& a, const std::vector<double>& b);
// other supplemental functions
double max(const double& x, const double& y);
int max(const int& x, const int& y);
long max(const long& x, const long& y);
#endif /* nonparametricSQG_hpp */
