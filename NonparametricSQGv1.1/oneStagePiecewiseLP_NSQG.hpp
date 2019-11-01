//
//  oneStagePiecewiseLP_NSQG.hpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/24/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#ifndef oneStagePiecewiseLP_NSQG_hpp
#define oneStagePiecewiseLP_NSQG_hpp

#include <stdio.h>
#include "nonparametricSQG.hpp"
#include "ioNonparametricModel.hpp"
#include "ioStochastic.hpp"
#include "ioNonparametricDB.hpp"

// One Stage Piecewise Linear Programming Solved by Nonparametric SQG
std::vector<double> oneStagePiecewiseLP_outputResults(const std::string& folder_path, int max_iterations, double initial_stepsize, std::vector<double> x_init);
// One Stage Piecewise Linear Programming Solved by Robust Nonparametric SQG
std::vector<double> robustOneStagePiecewiseLP_outputResults(const std::string& folder_path, int max_iterations, std::vector<double> x_init, int m, double M, double lowerBound, double upperBound);
// for debugging
std::vector<double> robustOneStagePiecewiseLP_outputResultsDEBUG(const std::string& folder_path, int max_iterations, std::vector<double> x_init, int m, double M, double lowerBound, double upperBound);
// Portforlio Optimization Validation by using true outcome (no prediction)
double oneStagePortfolioOptimizationValidation(const std::string& folder_path, const std::vector<double>& x_est, double W_trueResponse[], double cVaRRisk_alpha, double exchangeRate_lambda, double initialValue, int num_investments);
// One Stage Piecewise Linear Programming solved by SAA
std::vector<double> oneStagePiecewiseLP_SAA(const std::string& folder_path);
validationResult oneStagePiecewiseLP_SAA_validation(const std::string& folder_path, const std::vector<double>& x_candidate);
#endif /* oneStagePiecewiseLP_NSQG_hpp */
