//
//  twoStageLP_NSQG.hpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/23/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#ifndef twoStageLP_NSQG_hpp
#define twoStageLP_NSQG_hpp

#include <stdio.h>
#include <ilcplex/ilocplex.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <ctime>
#include "NonparametricSQG_dataType.hpp"
#include "nonparametricSQG.hpp"
#include "ioNonparametricDB.hpp"
#include "ioStochastic.hpp"
#include "ioNonparametricModel.hpp"

// random quantities only appear on the right hand side of constraints in the second stage problem
std::vector<double> twoStageLP_random_b(const std::string& folder_path, int max_iterations, double initial_stepsize, std::vector<double> x_init);
std::vector<double> twoStageLP_random_b_outputResults(const std::string& folder_path, int max_iterations, double initial_stepsize, std::vector<double> x_init);
// Robust Nonparametric SQG
std::vector<double> robustTwoStageLP_random_b_outputResults(const std::string& folder_path, int max_iterations, std::vector<double> x_init, double Dx, int m, double M);
// estimate cost via validation set
double twoStageLP_validation(const std::string& folder_path, const std::vector<double>& x_candidate);
double twoStageLP_validation_outputResults(const std::string& folder_path, const std::vector<double>& x_candidate);
double twoStageLP_secondStageCost(const std::vector<double>& x, const std::vector<double>& d, const std::vector<std::vector<double>>& De, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be, const std::vector<std::vector<double>>& Di, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi);
double twoStageLP_secondStageCost_inequality(const std::vector<double>& x, const std::vector<double>& d, const std::vector<std::vector<double>>& Di, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi);
double twoStageLP_secondStageCost_equality(const std::vector<double>& x, const std::vector<double>& d, const std::vector<std::vector<double>>& De, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be);
// estimate the value of nonparametric estimator
double twoStageLP_estimator(const std::string& folder_path, const std::vector<double>& x_candidate, int dataset_index);
// estimate the quality of candiate solution (v2) inclduing confidence interval
validationResult twoStageLP_validation_outputResultsV2(const std::string& folder_path, const std::vector<double>& x_candidate);
// find Dx
double Dx_estimateRectangle(const std::vector<double>&x, double lowerBound, double upperBound);
// test functions for debugging
void testTwoStageLPv1();
void testTwoStageLPValidation();
void projectionTest();
// iterative SAA
std::vector<double> twoStageLP_SAA(const std::string& folder_path);
#endif /* twoStageLP_NSQG_hpp */
