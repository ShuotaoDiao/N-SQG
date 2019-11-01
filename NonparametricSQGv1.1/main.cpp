//
//  main.cpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/23/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#include <iostream>
#include "twoStageLP_NSQG.hpp"
#include "oneStagePiecewiseLP_NSQG.hpp"
#include "testPlayground.hpp"

// generalized function for solving two stage shipment problem via NSQG
void NSQG_tssV2(int caseNumber, std::string method, std::vector<double> x, int maxIterations, double initStepsize) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/case" + std::to_string(caseNumber) + "/" + method;
    std::vector<double> x_candidate = twoStageLP_random_b_outputResults(folder_path, maxIterations, initStepsize, x);
    std::string ValidationFolder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment/validation";
    validationResult result = twoStageLP_validation_outputResultsV2(ValidationFolder_path, x_candidate);
    std::cout << "Average validation cost: " << result.mean << std::endl;
    std::cout << "Variance: " << result.variance << std::endl;
    std::cout << result.alpha << "% confidence interval: [" << result.CI_lower << ", " << result.CI_upper << "]" << std::endl;
}

void NSQG_tssV2_2(int caseNumber, std::string method, std::vector<double> x, int maxIterations, double initStepsize, int validationCaseNumber, int kNNValidationCaseNumber, int experimentNumber) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment" + std::to_string(experimentNumber) + "/replication" + std::to_string(caseNumber) + "/" + method;
    std::vector<double> x_candidate = twoStageLP_random_b_outputResults(folder_path, maxIterations, initStepsize, x);
    // Quality Estimation
    std::string validationFolder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment" + std::to_string(experimentNumber) +"/validation/validation" + std::to_string(validationCaseNumber);
    std::string kNNValidationFolder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment"+ std::to_string(experimentNumber) + "/kNNValidation" + std::to_string(kNNValidationCaseNumber);
    // estimate the quality of solution based on the true validation set
    validationResult trueResult = twoStageLP_validation_outputResultsV2(validationFolder_path, x_candidate);
    // estimate the quality of solution based on the kNN validation set
    validationResult kNNValidationResult = twoStageLP_validation_outputResultsV2(kNNValidationFolder_path, x_candidate);
    std::string validation_outputPath = folder_path + "/validationResults.txt";
    const char* validation_outputPath_const = validation_outputPath.c_str();
    std::fstream writeFile;
    writeFile.open(validation_outputPath_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    writeFile << "Solving Two Stage Shipment Problem by SQG-" + method + "\n";
    std::cout << "Solving Two Stage Shipment Problem by SQG-" + method + "\n";
    writeFile << "Initial Point        : ";
    std::cout << "Initial Point        : ";
    long x_size = x_candidate.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x[index] << ", ";
        std::cout << x[index] << ", ";
    }
    writeFile << x[x_size - 1] << "\n";
    std::cout << x[x_size - 1] << "\n";
    writeFile << "Initial Step Size    : " << initStepsize << "\n";
    std::cout << "Initial Step Size    : " << initStepsize << "\n";
    writeFile << "Number of Iterations : " << maxIterations << "\n";
    std::cout << "Number of Iterations : " << maxIterations << "\n";
    writeFile << "Estimating the quality of candidate solution\n";
    std::cout << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    std::cout << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution   : ";
    std::cout << "Candidate solution   : ";
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_candidate[index] << ", ";
        std::cout << x_candidate[index] << ", ";
    }
    writeFile << x_candidate[x_size - 1] << "\n";
    std::cout << x_candidate[x_size - 1] << "\n";
    writeFile << "========================True Validation Set========================\n";
    std::cout << "========================True Validation Set========================\n";
    writeFile << "Number of data points: " << trueResult.num_dataPoint << std::endl;
    std::cout << "Number of data points: " << trueResult.num_dataPoint << std::endl;
    writeFile << "Average cost         : " << trueResult.mean << std::endl;
    std::cout << "Average cost         : " << trueResult.mean << std::endl;
    writeFile << "Variance             : " << trueResult.variance << std::endl;
    std::cout << "Variance             : " << trueResult.variance << std::endl;
    writeFile << "Error in estimating the expected cost: " << sqrt(trueResult.variance / ((double)trueResult.num_dataPoint)) << std::endl;
    std::cout << "Error in estimating the expected cost: " << sqrt(trueResult.variance / ((double)trueResult.num_dataPoint)) << std::endl;
    writeFile << trueResult.alpha << "% CI Lower Bound: " << trueResult.CI_lower << std::endl;
    std::cout << trueResult.alpha << "% CI Lower Bound: " << trueResult.CI_lower << std::endl;
    writeFile << trueResult.alpha << "% CI Upper Bound: " << trueResult.CI_upper << std::endl;
    std::cout << trueResult.alpha << "% CI Upper Bound: " << trueResult.CI_upper << std::endl;
    writeFile << "========================kNN Validation Set=========================\n";
    std::cout << "========================kNN Validation Set=========================\n";
    writeFile << "Number of data points: " << kNNValidationResult.num_dataPoint << std::endl;
    std::cout << "Number of data points: " << kNNValidationResult.num_dataPoint << std::endl;
    writeFile << "Average cost         : " << kNNValidationResult.mean << std::endl;
    std::cout << "Average cost         : " << kNNValidationResult.mean << std::endl;
    writeFile << "Variance             : " << kNNValidationResult.variance << std::endl;
    std::cout << "Variance             : " << kNNValidationResult.variance << std::endl;
    writeFile << "Error in estimating the expected cost: " << sqrt(kNNValidationResult.variance / ((double)kNNValidationResult.num_dataPoint)) << std::endl;
    std::cout << "Error in estimating the expected cost: " << sqrt(kNNValidationResult.variance / ((double)kNNValidationResult.num_dataPoint)) << std::endl;
    writeFile << trueResult.alpha << "% CI Lower Bound: " << kNNValidationResult.CI_lower << std::endl;
    std::cout << trueResult.alpha << "% CI Lower Bound: " << kNNValidationResult.CI_lower << std::endl;
    writeFile << trueResult.alpha << "% CI Upper Bound: " << kNNValidationResult.CI_upper << std::endl;
    std::cout << trueResult.alpha << "% CI Upper Bound: " << kNNValidationResult.CI_upper << std::endl;
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    writeFile.close();
}

void SS_tssV2_2(int caseNumber, std::vector<double> x, int maxIterations, double initStepsize, int validationCaseNumber, int kNNValidationCaseNumber, int experimentNumber) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment" + std::to_string(experimentNumber) + "/replication" + std::to_string(caseNumber) + "/StochasticSubgradient";
    std::vector<double> x_candidate = twoStageLP_random_b_outputResults(folder_path, maxIterations, initStepsize, x);
    // Quality Estimation
    std::string validationFolder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment" + std::to_string(experimentNumber) + "/validation/validation" + std::to_string(validationCaseNumber);
    std::string kNNValidationFolder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment" + std::to_string(experimentNumber) + "/kNNValidation" + std::to_string(kNNValidationCaseNumber);
    // estimate the quality of solution based on the true validation set
    validationResult trueResult = twoStageLP_validation_outputResultsV2(validationFolder_path, x_candidate);
    // estimate the quality of solution based on the kNN validation set
    validationResult kNNValidationResult = twoStageLP_validation_outputResultsV2(kNNValidationFolder_path, x_candidate);
    std::string validation_outputPath = folder_path + "/validationResults.txt";
    const char* validation_outputPath_const = validation_outputPath.c_str();
    std::fstream writeFile;
    writeFile.open(validation_outputPath_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    writeFile << "Solving Two Stage Shipment Problem by Stochastic Subgradient Method\n";
    std::cout << "Solving Two Stage Shipment Problem by Stochastic Subgradient Method\n";
    writeFile << "Initial Point        : ";
    std::cout << "Initial Point        : ";
    long x_size = x_candidate.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x[index] << ", ";
        std::cout << x[index] << ", ";
    }
    writeFile << x[x_size - 1] << "\n";
    std::cout << x[x_size - 1] << "\n";
    writeFile << "Initial Step Size    : " << initStepsize << "\n";
    std::cout << "Initial Step Size    : " << initStepsize << "\n";
    writeFile << "Number of Iterations : " << maxIterations << "\n";
    std::cout << "Number of Iterations : " << maxIterations << "\n";
    writeFile << "Estimating the quality of candidate solution\n";
    std::cout << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    std::cout << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution   : ";
    std::cout << "Candidate solution   : ";
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_candidate[index] << ", ";
        std::cout << x_candidate[index] << ", ";
    }
    writeFile << x_candidate[x_size - 1] << "\n";
    std::cout << x_candidate[x_size - 1] << "\n";
    writeFile << "========================True Validation Set========================\n";
    std::cout << "========================True Validation Set========================\n";
    writeFile << "Number of data points: " << trueResult.num_dataPoint << std::endl;
    std::cout << "Number of data points: " << trueResult.num_dataPoint << std::endl;
    writeFile << "Average cost         : " << trueResult.mean << std::endl;
    std::cout << "Average cost         : " << trueResult.mean << std::endl;
    writeFile << "Variance             : " << trueResult.variance << std::endl;
    std::cout << "Variance             : " << trueResult.variance << std::endl;
    writeFile << "Error in estimating the expected cost: " << sqrt(trueResult.variance / ((double)trueResult.num_dataPoint)) << std::endl;
    std::cout << "Error in estimating the expected cost: " << sqrt(trueResult.variance / ((double)trueResult.num_dataPoint)) << std::endl;
    writeFile << trueResult.alpha << "% CI Lower Bound: " << trueResult.CI_lower << std::endl;
    std::cout << trueResult.alpha << "% CI Lower Bound: " << trueResult.CI_lower << std::endl;
    writeFile << trueResult.alpha << "% CI Upper Bound: " << trueResult.CI_upper << std::endl;
    std::cout << trueResult.alpha << "% CI Upper Bound: " << trueResult.CI_upper << std::endl;
    writeFile << "========================kNN Validation Set=========================\n";
    std::cout << "========================kNN Validation Set=========================\n";
    writeFile << "Number of data points: " << kNNValidationResult.num_dataPoint << std::endl;
    std::cout << "Number of data points: " << kNNValidationResult.num_dataPoint << std::endl;
    writeFile << "Average cost         : " << kNNValidationResult.mean << std::endl;
    std::cout << "Average cost         : " << kNNValidationResult.mean << std::endl;
    writeFile << "Variance             : " << kNNValidationResult.variance << std::endl;
    std::cout << "Variance             : " << kNNValidationResult.variance << std::endl;
    writeFile << "Error in estimating the expected cost: " << sqrt(kNNValidationResult.variance / ((double)kNNValidationResult.num_dataPoint)) << std::endl;
    std::cout << "Error in estimating the expected cost: " << sqrt(kNNValidationResult.variance / ((double)kNNValidationResult.num_dataPoint)) << std::endl;
    writeFile << trueResult.alpha << "% CI Lower Bound: " << kNNValidationResult.CI_lower << std::endl;
    std::cout << trueResult.alpha << "% CI Lower Bound: " << kNNValidationResult.CI_lower << std::endl;
    writeFile << trueResult.alpha << "% CI Upper Bound: " << kNNValidationResult.CI_upper << std::endl;
    std::cout << trueResult.alpha << "% CI Upper Bound: " << kNNValidationResult.CI_upper << std::endl;
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    writeFile.close();
}


// generalized function for solving two stage shipment problem via robust NSQG
void robustNSQG_tssV2(int caseNumber, std::string method, std::vector<double> x, int maxIterations) {
    double Dx = Dx_estimateRectangle(x, 0, 20);
    int m = 5;
    double M = 30;
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/case" + std::to_string(caseNumber) + "/" + method;
    std::vector<double> x_candidate = robustTwoStageLP_random_b_outputResults(folder_path, maxIterations, x, Dx, m, M);
    std::string ValidationFolder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/case" + std::to_string(caseNumber) + "/validation";
    validationResult result = twoStageLP_validation_outputResultsV2(ValidationFolder_path, x_candidate);
    std::cout << "Average validation cost: " << result.mean << std::endl;
    std::cout << "Variance: " << result.variance << std::endl;
    std::cout << result.alpha << "% confidence interval of expected cost: [" << result.CI_lower << ", " << result.CI_upper << "]" << std::endl;
}

void robustNSQG_tssV2_2(int caseNumber, std::string method, std::vector<double> x, int maxIterations, int validationCaseNumber, int kNNValidationCaseNumber, int experimentNumber) {
    double Dx = Dx_estimateRectangle(x, 0, 20);
    int m = 5;
    double M = 30;
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment" + std::to_string(experimentNumber) + "/replication" + std::to_string(caseNumber) + "/" + method;
    std::vector<double> x_candidate = robustTwoStageLP_random_b_outputResults(folder_path, maxIterations, x, Dx, m, M);
    //std::vector<double> x_candidate = x;
    // Quality Estimation
    std::string validationFolder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment" + std::to_string(experimentNumber) + "/validation/validation" + std::to_string(validationCaseNumber);
    std::string kNNValidationFolder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment" + std::to_string(experimentNumber) + "/kNNValidation" + std::to_string(kNNValidationCaseNumber);
    // estimate the quality of solution based on the true validation set
    validationResult trueResult = twoStageLP_validation_outputResultsV2(validationFolder_path, x_candidate);
    // estimate the quality of solution based on the kNN validation set
    validationResult kNNValidationResult = twoStageLP_validation_outputResultsV2(kNNValidationFolder_path, x_candidate);
    std::string validation_outputPath = folder_path + "/validationResults.txt";
    const char* validation_outputPath_const = validation_outputPath.c_str();
    std::fstream writeFile;
    writeFile.open(validation_outputPath_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    writeFile << "Solving Two Stage Shipment Problem by Robust SQG-" + method + "\n";
    std::cout << "Solving Two Stage Shipment Problem by Robust SQG-" + method + "\n";
    writeFile << "Initial Point        : ";
    std::cout << "Initial Point        : ";
    long x_size = x_candidate.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x[index] << ", ";
        std::cout << x[index] << ", ";
    }
    writeFile << x[x_size - 1] << "\n";
    std::cout << x[x_size - 1] << "\n";
    writeFile << "Number of Iterations (Outer loop): " << maxIterations << "\n";
    std::cout << "Number of Iterations (Outer loop): " << maxIterations << "\n";
    writeFile << "Dx                   : " << Dx << "\n";
    std::cout << "Dx                   : " << Dx << "\n";
    writeFile << "M                    : " << M << "\n";
    std::cout << "M                    : " << M << "\n";
    writeFile << "m                    : " << m << "\n";
    std::cout << "m                    : " << m << "\n";
    writeFile << "Estimating the quality of candidate solution\n";
    std::cout << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    std::cout << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution   : ";
    std::cout << "Candidate solution   : ";
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_candidate[index] << ", ";
        std::cout << x_candidate[index] << ", ";
    }
    writeFile << x_candidate[x_size - 1] << "\n";
    std::cout << x_candidate[x_size - 1] << "\n";
    writeFile << "========================True Validation Set========================\n";
    std::cout << "========================True Validation Set========================\n";
    writeFile << "Number of data points: " << trueResult.num_dataPoint << std::endl;
    std::cout << "Number of data points: " << trueResult.num_dataPoint << std::endl;
    writeFile << "Average cost         : " << trueResult.mean << std::endl;
    std::cout << "Average cost         : " << trueResult.mean << std::endl;
    writeFile << "Variance             : " << trueResult.variance << std::endl;
    std::cout << "Variance             : " << trueResult.variance << std::endl;
    writeFile << "Error in estimating the expected cost: " << sqrt(trueResult.variance / ((double)trueResult.num_dataPoint)) << std::endl;
    std::cout << "Error in estimating the expected cost: " << sqrt(trueResult.variance / ((double)trueResult.num_dataPoint)) << std::endl;
    writeFile << trueResult.alpha << "% CI Lower Bound: " << trueResult.CI_lower << std::endl;
    std::cout << trueResult.alpha << "% CI Lower Bound: " << trueResult.CI_lower << std::endl;
    writeFile << trueResult.alpha << "% CI Upper Bound: " << trueResult.CI_upper << std::endl;
    std::cout << trueResult.alpha << "% CI Upper Bound: " << trueResult.CI_upper << std::endl;
    writeFile << "========================kNN Validation Set=========================\n";
    std::cout << "========================kNN Validation Set=========================\n";
    writeFile << "Number of data points: " << kNNValidationResult.num_dataPoint << std::endl;
    std::cout << "Number of data points: " << kNNValidationResult.num_dataPoint << std::endl;
    writeFile << "Average cost         : " << kNNValidationResult.mean << std::endl;
    std::cout << "Average cost         : " << kNNValidationResult.mean << std::endl;
    writeFile << "Variance             : " << kNNValidationResult.variance << std::endl;
    std::cout << "Variance             : " << kNNValidationResult.variance << std::endl;
    writeFile << "Error in estimating the expected cost: " << sqrt(kNNValidationResult.variance / ((double)kNNValidationResult.num_dataPoint)) << std::endl;
    std::cout << "Error in estimating the expected cost: " << sqrt(kNNValidationResult.variance / ((double)kNNValidationResult.num_dataPoint)) << std::endl;
    writeFile << trueResult.alpha << "% CI Lower Bound: " << kNNValidationResult.CI_lower << std::endl;
    std::cout << trueResult.alpha << "% CI Lower Bound: " << kNNValidationResult.CI_lower << std::endl;
    writeFile << trueResult.alpha << "% CI Upper Bound: " << kNNValidationResult.CI_upper << std::endl;
    std::cout << trueResult.alpha << "% CI Upper Bound: " << kNNValidationResult.CI_upper << std::endl;
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    writeFile.close();
}

void robustSS_tssV2_2(int caseNumber, std::vector<double> x, int maxIterations, int validationCaseNumber, int kNNValidationCaseNumber, int experimentNumber) {
    double Dx = Dx_estimateRectangle(x, 0, 20);
    int m = 5;
    double M = 30;
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment" + std::to_string(experimentNumber) + "/replication" + std::to_string(caseNumber) + "/StochasticSubgradient";
    std::vector<double> x_candidate = robustTwoStageLP_random_b_outputResults(folder_path, maxIterations, x, Dx, m, M);
    // Quality Estimation
    std::string validationFolder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment3/validation/validation" + std::to_string(validationCaseNumber);
    std::string kNNValidationFolder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment3/kNNValidation" + std::to_string(kNNValidationCaseNumber);
    // estimate the quality of solution based on the true validation set
    validationResult trueResult = twoStageLP_validation_outputResultsV2(validationFolder_path, x_candidate);
    // estimate the quality of solution based on the kNN validation set
    validationResult kNNValidationResult = twoStageLP_validation_outputResultsV2(kNNValidationFolder_path, x_candidate);
    std::string validation_outputPath = folder_path + "/validationResults.txt";
    const char* validation_outputPath_const = validation_outputPath.c_str();
    std::fstream writeFile;
    writeFile.open(validation_outputPath_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    writeFile << "Solving Two Stage Shipment Problem by Robust Stochastic Subgradient Method\n";
    std::cout << "Solving Two Stage Shipment Problem by Robust Stochastic Subgradient Method\n";
    writeFile << "Initial Point        : ";
    std::cout << "Initial Point        : ";
    long x_size = x_candidate.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x[index] << ", ";
        std::cout << x[index] << ", ";
    }
    writeFile << x[x_size - 1] << "\n";
    std::cout << x[x_size - 1] << "\n";
    writeFile << "Number of Iterations (Outer loop): " << maxIterations << "\n";
    std::cout << "Number of Iterations (Outer loop): " << maxIterations << "\n";
    writeFile << "Dx                   : " << Dx << "\n";
    std::cout << "Dx                   : " << Dx << "\n";
    writeFile << "M                    : " << M << "\n";
    std::cout << "M                    : " << M << "\n";
    writeFile << "m                    : " << m << "\n";
    std::cout << "m                    : " << m << "\n";
    writeFile << "Estimating the quality of candidate solution\n";
    std::cout << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    std::cout << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution   : ";
    std::cout << "Candidate solution   : ";
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_candidate[index] << ", ";
        std::cout << x_candidate[index] << ", ";
    }
    writeFile << x_candidate[x_size - 1] << "\n";
    std::cout << x_candidate[x_size - 1] << "\n";
    writeFile << "========================True Validation Set========================\n";
    std::cout << "========================True Validation Set========================\n";
    writeFile << "Number of data points: " << trueResult.num_dataPoint << std::endl;
    std::cout << "Number of data points: " << trueResult.num_dataPoint << std::endl;
    writeFile << "Average cost         : " << trueResult.mean << std::endl;
    std::cout << "Average cost         : " << trueResult.mean << std::endl;
    writeFile << "Variance             : " << trueResult.variance << std::endl;
    std::cout << "Variance             : " << trueResult.variance << std::endl;
    writeFile << "Error in estimating the expected cost: " << sqrt(trueResult.variance / ((double)trueResult.num_dataPoint)) << std::endl;
    std::cout << "Error in estimating the expected cost: " << sqrt(trueResult.variance / ((double)trueResult.num_dataPoint)) << std::endl;
    writeFile << trueResult.alpha << "% CI Lower Bound: " << trueResult.CI_lower << std::endl;
    std::cout << trueResult.alpha << "% CI Lower Bound: " << trueResult.CI_lower << std::endl;
    writeFile << trueResult.alpha << "% CI Upper Bound: " << trueResult.CI_upper << std::endl;
    std::cout << trueResult.alpha << "% CI Upper Bound: " << trueResult.CI_upper << std::endl;
    writeFile << "========================kNN Validation Set=========================\n";
    std::cout << "========================kNN Validation Set=========================\n";
    writeFile << "Number of data points: " << kNNValidationResult.num_dataPoint << std::endl;
    std::cout << "Number of data points: " << kNNValidationResult.num_dataPoint << std::endl;
    writeFile << "Average cost         : " << kNNValidationResult.mean << std::endl;
    std::cout << "Average cost         : " << kNNValidationResult.mean << std::endl;
    writeFile << "Variance             : " << kNNValidationResult.variance << std::endl;
    std::cout << "Variance             : " << kNNValidationResult.variance << std::endl;
    writeFile << "Error in estimating the expected cost: " << sqrt(kNNValidationResult.variance / ((double)kNNValidationResult.num_dataPoint)) << std::endl;
    std::cout << "Error in estimating the expected cost: " << sqrt(kNNValidationResult.variance / ((double)kNNValidationResult.num_dataPoint)) << std::endl;
    writeFile << trueResult.alpha << "% CI Lower Bound: " << kNNValidationResult.CI_lower << std::endl;
    std::cout << trueResult.alpha << "% CI Lower Bound: " << kNNValidationResult.CI_lower << std::endl;
    writeFile << trueResult.alpha << "% CI Upper Bound: " << kNNValidationResult.CI_upper << std::endl;
    std::cout << trueResult.alpha << "% CI Upper Bound: " << kNNValidationResult.CI_upper << std::endl;
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    writeFile.close();
}

void SAA_tsskNNValidation(int kNNValidationCaseNumber){
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment3/kNNValidation" + std::to_string(kNNValidationCaseNumber);
    std::vector<double> x_candidate = twoStageLP_SAA(folder_path);
}
//
void SAA_tssTrueValidation(int validationCaseNumber){
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment3/validation/validation" + std::to_string(validationCaseNumber);
    std::vector<double> x_candidate = twoStageLP_SAA(folder_path);// output results
}

// SAA for solving two stage shipment problem
void SAA_tss(int caseNumber, int num_datasets, int validationCaseNumber, int kNNValidationCaseNumber) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment3/replication" + std::to_string(caseNumber) + "/SAA_numDatasets" + std::to_string(num_datasets);
    std::vector<double> x_candidate = twoStageLP_SAA(folder_path);
    std::string validationFolder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment3/validation/validation" + std::to_string(validationCaseNumber);
    std::string kNNValidationFolder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment3/kNNValidation" + std::to_string(kNNValidationCaseNumber);
    // estimate the quality of solution based on the true validation set
    validationResult trueResult = twoStageLP_validation_outputResultsV2(validationFolder_path, x_candidate);
    // estimate the quality of solution based on the kNN validation set
    validationResult kNNValidationResult = twoStageLP_validation_outputResultsV2(kNNValidationFolder_path, x_candidate);
    std::string validation_outputPath = folder_path + "/validationResults.txt";
    const char* validation_outputPath_const = validation_outputPath.c_str();
    std::fstream writeFile;
    writeFile.open(validation_outputPath_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    writeFile << "Solving Two Stage Shipment Problem by SAA\n";
    std::cout << "Solving Two Stage Shipment Problem by SAA\n";
    writeFile << "Estimating the quality of candidate solution\n";
    std::cout << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    std::cout << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution   : ";
    std::cout << "Candidate solution   : ";
    long x_size = x_candidate.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_candidate[index] << ", ";
        std::cout << x_candidate[index] << ", ";
    }
    writeFile << x_candidate[x_size - 1] << "\n";
    std::cout << x_candidate[x_size - 1] << "\n";
    writeFile << "========================True Validation Set========================\n";
    std::cout << "========================True Validation Set========================\n";
    writeFile << "Number of data points: " << trueResult.num_dataPoint << std::endl;
    std::cout << "Number of data points: " << trueResult.num_dataPoint << std::endl;
    writeFile << "Average cost         : " << trueResult.mean << std::endl;
    std::cout << "Average cost         : " << trueResult.mean << std::endl;
    writeFile << "Variance             : " << trueResult.variance << std::endl;
    std::cout << "Variance             : " << trueResult.variance << std::endl;
    writeFile << "Error in estimating the expected cost: " << sqrt(trueResult.variance / ((double)trueResult.num_dataPoint)) << std::endl;
    std::cout << "Error in estimating the expected cost: " << sqrt(trueResult.variance / ((double)trueResult.num_dataPoint)) << std::endl;
    writeFile << trueResult.alpha << "% CI Lower Bound: " << trueResult.CI_lower << std::endl;
    std::cout << trueResult.alpha << "% CI Lower Bound: " << trueResult.CI_lower << std::endl;
    writeFile << trueResult.alpha << "% CI Upper Bound: " << trueResult.CI_upper << std::endl;
    std::cout << trueResult.alpha << "% CI Upper Bound: " << trueResult.CI_upper << std::endl;
    writeFile << "========================kNN Validation Set=========================\n";
    std::cout << "========================kNN Validation Set=========================\n";
    writeFile << "Number of data points: " << kNNValidationResult.num_dataPoint << std::endl;
    std::cout << "Number of data points: " << kNNValidationResult.num_dataPoint << std::endl;
    writeFile << "Average cost         : " << kNNValidationResult.mean << std::endl;
    std::cout << "Average cost         : " << kNNValidationResult.mean << std::endl;
    writeFile << "Variance             : " << kNNValidationResult.variance << std::endl;
    std::cout << "Variance             : " << kNNValidationResult.variance << std::endl;
    writeFile << "Error in estimating the expected cost: " << sqrt(kNNValidationResult.variance / ((double)kNNValidationResult.num_dataPoint)) << std::endl;
    std::cout << "Error in estimating the expected cost: " << sqrt(kNNValidationResult.variance / ((double)kNNValidationResult.num_dataPoint)) << std::endl;
    writeFile << trueResult.alpha << "% CI Lower Bound: " << kNNValidationResult.CI_lower << std::endl;
    std::cout << trueResult.alpha << "% CI Lower Bound: " << kNNValidationResult.CI_lower << std::endl;
    writeFile << trueResult.alpha << "% CI Upper Bound: " << kNNValidationResult.CI_upper << std::endl;
    std::cout << trueResult.alpha << "% CI Upper Bound: " << kNNValidationResult.CI_upper << std::endl;
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    writeFile.close();
}

// calculate the error term
void tss_errorCalculation(int caseNumber, std::string method, std::vector<double> x, int experimentNumber, int dataset_index) {
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment" + std::to_string(experimentNumber) + "/replication" + std::to_string(caseNumber) + "/" + method;
    double estimatedCost = twoStageLP_estimator(folder_path, x, dataset_index);
    std::cout << "**********************************************************************\n";
    std::cout << "Folder Path                              : " << folder_path << std::endl;
    std::cout << "Given Point                              : ";
    for (int x_index = 0; x_index < x.size() - 1; ++x_index) {
        std::cout << x[x_index] << ", ";
    }
    std::cout << x[x.size() - 1] << std::endl;
    std::cout << "Cost of Estimated Function at Given Point: " << estimatedCost << std::endl;
    std::cout << "**********************************************************************\n";
}

void tss_errorCalculationV2(int experimentNumber, std::string method, std::vector<double> x, int dataset_index) {
    validationResult result;
    double mean = 0;
    double tempV = 0; // for calculating variance
    for (int case_index = 1; case_index < 31; ++case_index) {
        std::string folder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment" + std::to_string(experimentNumber) + "/replication" + std::to_string(case_index) + "/" + method;
        double estimatedCost = twoStageLP_estimator(folder_path, x, dataset_index);
        mean += estimatedCost / 30.0;
        tempV = tempV * (case_index - 1) / case_index + estimatedCost * estimatedCost / case_index;
        std::cout << "**********************************************************************\n";
        std::cout << "Folder Path                              : " << folder_path << std::endl;
        std::cout << "Given Point                              : ";
        for (int x_index = 0; x_index < x.size() - 1; ++x_index) {
            std::cout << x[x_index] << ", ";
        }
        std::cout << x[x.size() - 1] << std::endl;
        std::cout << "Cost of Estimated Function at Given Point: " << estimatedCost << std::endl;
        std::cout << "**********************************************************************\n";
    }
    result.mean = mean;
    result.variance = tempV * 30.0 / 29.0 - mean * mean * 30.0 / 29.0;
    double error = sqrt(result.variance / 30.0);
    double halfMargin = error * result.Zalpha;
    result.CI_lower = mean - halfMargin;
    result.CI_upper = mean + halfMargin;
    std::string outputPath = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment" + std::to_string(experimentNumber) + "/errorEstimation.txt";
    const char* outputPath_const = outputPath.c_str();
    std::fstream writeFile;
    writeFile.open(outputPath_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    writeFile << "Candidate Point      :";
    for (int x_index = 0; x_index < x.size() - 1; ++x_index) {
        writeFile << x[x_index] << ", ";
    }
    writeFile << x[x.size() - 1] << std::endl;
    writeFile << "Estimating the bias of " << method << " estimator" << std::endl;
    std::cout << "Estimating the bias of " << method << " estimator" << std::endl;
    writeFile << "Dataset Index        : " << dataset_index << std::endl;
    std::cout << "Dataset Index        : " << dataset_index << std::endl;
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    std::cout << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Average cost         : " << result.mean << std::endl;
    std::cout << "Average cost         : " << result.mean << std::endl;
    writeFile << "Variance             : " << result.variance << std::endl;
    std::cout << "Variance             : " << result.variance << std::endl;
    writeFile << "Error in estimating the expected cost: " << error << std::endl;
    std::cout << "Error in estimating the expected cost: " << error << std::endl;
    writeFile << result.alpha << "% CI Lower Bound: " << result.CI_lower << std::endl;
    std::cout << result.alpha << "% CI Lower Bound: " << result.CI_lower << std::endl;
    writeFile << result.alpha << "% CI Upper Bound: " << result.CI_upper << std::endl;
    std::cout << result.alpha << "% CI Upper Bound: " << result.CI_upper << std::endl;
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
}

void tss_errorCalculationV3(int experimentNumber, std::string method, std::vector<double> x, int dataset_index, int num_cases) {
    validationResult result;
    double mean = 0;
    double tempV = 0; // for calculating variance
    for (int case_index = 1; case_index < num_cases + 1; ++case_index) {
        std::string folder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment" + std::to_string(experimentNumber) + "/replication" + std::to_string(case_index) + "/" + method;
        double estimatedCost = twoStageLP_estimator(folder_path, x, dataset_index);
        mean += estimatedCost / ((double) num_cases);
        tempV = tempV * (case_index - 1) / case_index + estimatedCost * estimatedCost / case_index;
        std::cout << "**********************************************************************\n";
        std::cout << "Folder Path                              : " << folder_path << std::endl;
        std::cout << "Given Point                              : ";
        for (int x_index = 0; x_index < x.size() - 1; ++x_index) {
            std::cout << x[x_index] << ", ";
        }
        std::cout << x[x.size() - 1] << std::endl;
        std::cout << "Cost of Estimated Function at Given Point: " << estimatedCost << std::endl;
        std::cout << "**********************************************************************\n";
    }
    result.mean = mean;
    result.variance = tempV * ((double) num_cases) / ((double) num_cases - 1) - mean * mean * ((double) num_cases) / ((double) num_cases - 1);
    double error = sqrt(result.variance / ((double) num_cases));
    double halfMargin = error * result.Zalpha;
    result.CI_lower = mean - halfMargin;
    result.CI_upper = mean + halfMargin;
    std::string outputPath = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/experiment" + std::to_string(experimentNumber) + "/errorEstimationV3.txt";
    const char* outputPath_const = outputPath.c_str();
    std::fstream writeFile;
    writeFile.open(outputPath_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    std::cout << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Number of replications: " << num_cases << std::endl;
    std::cout << "Number of replications: " << num_cases << std::endl;
    writeFile << "Candidate Point      :";
    for (int x_index = 0; x_index < x.size() - 1; ++x_index) {
        writeFile << x[x_index] << ", ";
    }
    writeFile << x[x.size() - 1] << std::endl;
    writeFile << "Estimating the bias of " << method << " estimator" << std::endl;
    std::cout << "Estimating the bias of " << method << " estimator" << std::endl;
    writeFile << "Dataset Index        : " << dataset_index << std::endl;
    std::cout << "Dataset Index        : " << dataset_index << std::endl;
    writeFile << "Average cost         : " << result.mean << std::endl;
    std::cout << "Average cost         : " << result.mean << std::endl;
    writeFile << "Variance             : " << result.variance << std::endl;
    std::cout << "Variance             : " << result.variance << std::endl;
    writeFile << "Error in estimating the expected cost: " << error << std::endl;
    std::cout << "Error in estimating the expected cost: " << error << std::endl;
    writeFile << result.alpha << "% CI Lower Bound: " << result.CI_lower << std::endl;
    std::cout << result.alpha << "% CI Lower Bound: " << result.CI_lower << std::endl;
    writeFile << result.alpha << "% CI Upper Bound: " << result.CI_upper << std::endl;
    std::cout << result.alpha << "% CI Upper Bound: " << result.CI_upper << std::endl;
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    writeFile.close();
}


// =================================================================================================
// one stage piecewise linear portfolio optimization problem via NSQG
void NSQG_portfolioOptimization(std::string method, std::vector<double> x, int maxIterations, double initialStepsize) {
    //std::string caseNumber = "2";
    std::string targetTime = "2015";
    // true response
    double W_trueResponse[] = {1, 139.75705, 47.207386, 119.690994, 141.210007, 235.925293, 59.700073, 85.5019, 164.313843};
    //std::string folder_path = "/Users/sonny/Documents/S&P100/case" + caseNumber + "/" + method;
    std::string folder_path = "/Users/sonny/Documents/S&P100/8stocks/targetTime_" + targetTime + "/" + method;
    std::string validationFolder_path = "/Users/sonny/Documents/S&P100/8stocks/targetTime_" + targetTime + "/SAA_kNN";
    std::vector<double> x_candidate = oneStagePiecewiseLP_outputResults(folder_path, maxIterations, initialStepsize, x);
    std::cout << "Candidate Solution for Portfolio Optimization by Nonparametric SQG" << std::endl;
    for (int x_index = 0; x_index < x_candidate.size(); ++x_index) {
        std::cout << "x[" << x_index << "] = " << x_candidate[x_index] << std::endl;
    }
    int num_investments = (sizeof(W_trueResponse) / sizeof(W_trueResponse[0]));
    std::cout << num_investments << std::endl;
    double initialValue = 100000;
    double cVaRRisk_alpha = 0.1;
    double exchangeRate_lambda = 0.2;
    std::cout << "Start Measuring the Quality of the Solution" << std::endl;
    // get the lowest cost of target day
    double lowestTrueCost = oneStagePortfolioOptimizationValidation(folder_path, x_candidate, W_trueResponse, cVaRRisk_alpha, exchangeRate_lambda, initialValue, num_investments);
    // get the objective cost
    double trueCost = 0;
    trueCost = x_candidate[0];
    int base_index = 1;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        trueCost -= exchangeRate_lambda * W_trueResponse[investment_index] * x_candidate[base_index +investment_index];
    }
    // obtain the value in the piecewise linear function
    double PL_value = 0;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        PL_value += -W_trueResponse[investment_index] * x_candidate[base_index + investment_index];
    }
    PL_value += initialValue - x_candidate[0];
    if (PL_value > 0) {
        trueCost += (1.0 / (1.0 - cVaRRisk_alpha)) * PL_value;
    }
    validationResult result = oneStagePiecewiseLP_SAA_validation(validationFolder_path, x_candidate);
    // write the results of quality of candidate solution
    std::string outputResults_path = folder_path + "/validationResults.txt";
    const char* outputResults_path_const = outputResults_path.c_str();
    std::fstream writeFile;
    writeFile.open(outputResults_path_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "***************************************************\n";
    writeFile << "Nonparametric SQG-" + method + " has been used" << std::endl;
    writeFile << "Number of iterations: " << maxIterations << std::endl;
    writeFile << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution: ";
    long x_size = x_candidate.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_candidate[index] << ", ";
    }
    writeFile << x_candidate[x_size - 1] << "\n";
    writeFile << "True Cost (by using the true outcome): " << trueCost << std::endl;
    std::cout << "True Cost (by using the true outcome): " << trueCost << std::endl;
    writeFile << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    std::cout << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "Validation Set Generated by Getting the k Nearest Neighbors of the whole Training Set" << std::endl;
    writeFile << "Mean of cost: " << result.mean << std::endl;
    writeFile << "Variance of cost: " << result.variance << std::endl;
    writeFile << result.alpha << "% CI Lower Bound: " << result.CI_lower << std::endl;
    writeFile << result.alpha << "% CI Upper Bound: " << result.CI_upper << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "***************************************************\n";
    writeFile.close();
}

// one stage piecewise linear portfolio optimization problem via NSQG
void NSQG_portfolioOptimizationV2(std::string method, std::vector<double> x, int maxIterations, double initialStepsize, std::string targetTime, double W_trueResponse[], int num_investments) {
    //std::string folder_path = "/Users/sonny/Documents/S&P100/case" + caseNumber + "/" + method;
    std::string folder_path = "/Users/sonny/Documents/S&P100/8stocks2.0/targetTime_" + targetTime + "/" + method;
    std::string validationFolder_path = "/Users/sonny/Documents/S&P100/8stocks2.0/targetTime_" + targetTime + "/SAA_kNN";
    std::vector<double> x_candidate = oneStagePiecewiseLP_outputResults(folder_path, maxIterations, initialStepsize, x);
    std::cout << "Candidate Solution for Portfolio Optimization by Nonparametric SQG" << std::endl;
    for (int x_index = 0; x_index < x_candidate.size(); ++x_index) {
        std::cout << "x[" << x_index << "] = " << x_candidate[x_index] << std::endl;
    }
    std::cout << num_investments << std::endl;
    double initialValue = 100000;
    double cVaRRisk_alpha = 0.1;
    double exchangeRate_lambda = 0.2;
    std::cout << "Start Measuring the Quality of the Solution" << std::endl;
    // get the lowest cost of target day
    double lowestTrueCost = oneStagePortfolioOptimizationValidation(folder_path, x_candidate, W_trueResponse, cVaRRisk_alpha, exchangeRate_lambda, initialValue, num_investments);
    // get the objective cost
    double trueCost = 0;
    trueCost = x_candidate[0];
    int base_index = 1;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        trueCost -= exchangeRate_lambda * W_trueResponse[investment_index] * x_candidate[base_index +investment_index];
    }
    // obtain the value in the piecewise linear function
    double PL_value = 0;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        PL_value += -W_trueResponse[investment_index] * x_candidate[base_index + investment_index];
    }
    PL_value += initialValue - x_candidate[0];
    if (PL_value > 0) {
        trueCost += (1.0 / (1.0 - cVaRRisk_alpha)) * PL_value;
    }
    validationResult result = oneStagePiecewiseLP_SAA_validation(validationFolder_path, x_candidate);
    // write the results of quality of candidate solution
    std::string outputResults_path = folder_path + "/validationResults.txt";
    const char* outputResults_path_const = outputResults_path.c_str();
    std::fstream writeFile;
    writeFile.open(outputResults_path_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "***************************************************\n";
    writeFile << "Nonparametric SQG-" + method + " has been used" << std::endl;
    writeFile << "Number of iterations: " << maxIterations << std::endl;
    writeFile << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution: ";
    long x_size = x_candidate.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_candidate[index] << ", ";
    }
    writeFile << x_candidate[x_size - 1] << "\n";
    writeFile << "True Cost (by using the true outcome): " << trueCost << std::endl;
    std::cout << "True Cost (by using the true outcome): " << trueCost << std::endl;
    writeFile << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    std::cout << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "Validation Set Generated by Getting the k Nearest Neighbors of the whole Training Set" << std::endl;
    writeFile << "Mean of cost: " << result.mean << std::endl;
    writeFile << "Variance of cost: " << result.variance << std::endl;
    writeFile << result.alpha << "% CI Lower Bound: " << result.CI_lower << std::endl;
    writeFile << result.alpha << "% CI Upper Bound: " << result.CI_upper << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "***************************************************\n";
    writeFile.close();
}

// one stage piecewise linear portfolio optimization problem via robust NSQG
void robustNSQG_portfolioOptimization(std::string method, std::vector<double> x, int maxIterations) {
    double M = 307.6386;
    double lowerBound = 0;
    double upperBound = 200;
    double Dx = Dx_estimateRectangle(x, 0, 200);
    //std::string caseNumber = "2";
    std::string targetTime = "2015";
    // true response
    double W_trueResponse[] = {1, 139.75705, 47.207386, 119.690994, 141.210007, 235.925293, 59.700073, 85.5019, 164.313843};
    std::cout << "Dx: " << Dx << std::endl;
    int m = 5;
    //std::string folder_path = "/Users/sonny/Documents/S&P100/case" + caseNumber + "/" + method;
    std::string folder_path = "/Users/sonny/Documents/S&P100/8stocks/targetTime_" + targetTime + "/" + method;
    std::string validationFolder_path = "/Users/sonny/Documents/S&P100/8stocks/targetTime_" + targetTime + "/SAA_kNN";
    // solve by using robust NSQG
    std::vector<double> x_candidate = robustOneStagePiecewiseLP_outputResults(folder_path, maxIterations, x, m, M, lowerBound, upperBound);
    int num_investments = (sizeof(W_trueResponse) / sizeof(W_trueResponse[0]));
    std::cout << num_investments << std::endl;
    double initialValue = 100000;
    double cVaRRisk_alpha = 0.1;
    double exchangeRate_lambda = 0.2;
    std::cout << "Start Measuring the Quality of the Solution" << std::endl;
    // get the lowest cost of target day
    double lowestTrueCost = oneStagePortfolioOptimizationValidation(folder_path, x_candidate, W_trueResponse, cVaRRisk_alpha, exchangeRate_lambda, initialValue, num_investments);
    // get the objective cost
    double trueCost = 0;
    trueCost = x_candidate[0];
    int base_index = 1;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        trueCost -= exchangeRate_lambda * W_trueResponse[investment_index] * x_candidate[base_index +investment_index];
    }
    // obtain the value in the piecewise linear function
    double PL_value = 0;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        PL_value += -W_trueResponse[investment_index] * x_candidate[base_index + investment_index];
    }
    PL_value += initialValue - x_candidate[0];
    if (PL_value > 0) {
        trueCost += (1.0 / (1.0 - cVaRRisk_alpha)) * PL_value;
    }
    validationResult result = oneStagePiecewiseLP_SAA_validation(validationFolder_path, x_candidate);
    // write the results of quality of candidate solution
    std::string outputResults_path = folder_path + "/validationResults.txt";
    const char* outputResults_path_const = outputResults_path.c_str();
    std::fstream writeFile;
    writeFile.open(outputResults_path_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "***************************************************\n";
    writeFile << "Robust Nonparametric SQG-" + method + " has been used" << std::endl;
    writeFile << "Number of outer loops: " << maxIterations << std::endl;
    writeFile << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution: ";
    long x_size = x_candidate.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_candidate[index] << ", ";
    }
    writeFile << x_candidate[x_size - 1] << "\n";
    writeFile << "True Cost (by using the true outcome): " << trueCost << std::endl;
    std::cout << "True Cost (by using the true outcome): " << trueCost << std::endl;
    writeFile << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    std::cout << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "Validation Set Generated by Getting the k Nearest Neighbors of the whole Training Set" << std::endl;
    writeFile << "Mean of cost: " << result.mean << std::endl;
    writeFile << "Variance of cost: " << result.variance << std::endl;
    writeFile << result.alpha << "% CI Lower Bound: " << result.CI_lower << std::endl;
    writeFile << result.alpha << "% CI Upper Bound: " << result.CI_upper << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "***************************************************\n";
    writeFile.close();
}

// one stage piecewise linear portfolio optimization problem via robust NSQG
void robustNSQG_portfolioOptimizationV2(std::string method, std::vector<double> x, int maxIterations, std::string targetTime, double W_trueResponse[], int num_investments) {
    double M = 307.6386;
    double lowerBound = 0;
    double upperBound = 200;
    double Dx = Dx_estimateRectangle(x, 0, 200);
    std::cout << "Dx: " << Dx << std::endl;
    int m = 5;
    //std::string folder_path = "/Users/sonny/Documents/S&P100/case" + caseNumber + "/" + method;
    std::string folder_path = "/Users/sonny/Documents/S&P100/8stocks2.0/targetTime_" + targetTime + "/" + method;
    std::string validationFolder_path = "/Users/sonny/Documents/S&P100/8stocks2.0/targetTime_" + targetTime + "/SAA_kNN";
    // solve by using robust NSQG
    std::vector<double> x_candidate = robustOneStagePiecewiseLP_outputResults(folder_path, maxIterations, x, m, M, lowerBound, upperBound);
    std::cout << num_investments << std::endl;
    double initialValue = 100000;
    double cVaRRisk_alpha = 0.1;
    double exchangeRate_lambda = 0.2;
    std::cout << "Start Measuring the Quality of the Solution" << std::endl;
    // get the lowest cost of target day
    double lowestTrueCost = oneStagePortfolioOptimizationValidation(folder_path, x_candidate, W_trueResponse, cVaRRisk_alpha, exchangeRate_lambda, initialValue, num_investments);
    // get the objective cost
    double trueCost = 0;
    trueCost = x_candidate[0];
    int base_index = 1;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        trueCost -= exchangeRate_lambda * W_trueResponse[investment_index] * x_candidate[base_index +investment_index];
    }
    // obtain the value in the piecewise linear function
    double PL_value = 0;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        PL_value += -W_trueResponse[investment_index] * x_candidate[base_index + investment_index];
    }
    PL_value += initialValue - x_candidate[0];
    if (PL_value > 0) {
        trueCost += (1.0 / (1.0 - cVaRRisk_alpha)) * PL_value;
    }
    validationResult result = oneStagePiecewiseLP_SAA_validation(validationFolder_path, x_candidate);
    // write the results of quality of candidate solution
    std::string outputResults_path = folder_path + "/validationResults.txt";
    const char* outputResults_path_const = outputResults_path.c_str();
    std::fstream writeFile;
    writeFile.open(outputResults_path_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "***************************************************\n";
    writeFile << "Robust Nonparametric SQG-" + method + " has been used" << std::endl;
    writeFile << "Number of outer loops: " << maxIterations << std::endl;
    writeFile << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution: ";
    long x_size = x_candidate.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_candidate[index] << ", ";
    }
    writeFile << x_candidate[x_size - 1] << "\n";
    writeFile << "True Cost (by using the true outcome): " << trueCost << std::endl;
    std::cout << "True Cost (by using the true outcome): " << trueCost << std::endl;
    writeFile << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    std::cout << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "Validation Set Generated by Getting the k Nearest Neighbors of the whole Training Set" << std::endl;
    writeFile << "Mean of cost: " << result.mean << std::endl;
    writeFile << "Variance of cost: " << result.variance << std::endl;
    writeFile << result.alpha << "% CI Lower Bound: " << result.CI_lower << std::endl;
    writeFile << result.alpha << "% CI Upper Bound: " << result.CI_upper << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "***************************************************\n";
    writeFile.close();
}


// debugging robust Nonparametric SQG
void robustNSQG_portfolioOptimizationDEBUG(std::string method, std::vector<double> x, int maxIterations) {
    double M = 307.6386;
    double lowerBound = 0;
    double upperBound = 200;
    double Dx = Dx_estimateRectangle(x, 0, 200);
    //std::string caseNumber = "2";
    std::string targetTime = "2014";
    // true response
    double W_trueResponse[] = {1, 139.184799, 46.92799, 119.575211, 138.410004, 235.625885, 60.185062, 85.328598, 162.937164};
    std::cout << "Dx: " << Dx << std::endl;
    int m = 5;
    //std::string folder_path = "/Users/sonny/Documents/S&P100/case" + caseNumber + "/" + method;
    std::string folder_path = "/Users/sonny/Documents/S&P100/8stocks/targetTime_" + targetTime + "/" + method;
    std::string validationFolder_path = "/Users/sonny/Documents/S&P100/8stocks/targetTime_" + targetTime + "/SAA_kNN";
    // solve by using robust NSQG
    std::vector<double> x_candidate = robustOneStagePiecewiseLP_outputResultsDEBUG(folder_path, maxIterations, x, m, M, lowerBound, upperBound);
}

void oneStageValidationTest() {
    std::vector<double> x_candidate(28,0);
    std::string method = "kNN";
    std::string targetTime = "2013";
    std::string folder_path = "/Users/sonny/Documents/S&P100/8stocks/targetTime_" + targetTime + "/" + method;
    double W_trueResponse[] = {1, 139.689148, 46.83165, 119.334007, 141.440002, 234.766342, 59.452835, 85.145683, 162.399719};
    int num_investments = (sizeof(W_trueResponse) / sizeof(W_trueResponse[0]));
    double initialValue = 100000;
    double cVaRRisk_alpha = 0.1;
    double exchangeRate_lambda = 0.2;
    // get the lowest cost of target day
    double lowestTrueCost = oneStagePortfolioOptimizationValidation(folder_path, x_candidate, W_trueResponse, cVaRRisk_alpha, exchangeRate_lambda, initialValue, num_investments);
    std::cout << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
}

void oneStagePortfolioOtimization_SAA() {
    std::string targetTime = "2015";
    double W_trueResponse[] = {1, 139.75705, 47.207386, 119.690994, 141.210007, 235.925293, 59.700073, 85.5019, 164.313843};
    int num_investments = (sizeof(W_trueResponse) / sizeof(W_trueResponse[0]));
    double initialValue = 100000;
    double cVaRRisk_alpha = 0.1;
    double exchangeRate_lambda = 0.2;
    std::string folder_path = "/Users/sonny/Documents/S&P100/8stocks/targetTime_" + targetTime + "/SAA";
    std::string outputResults_path = folder_path + "/validationResults.txt";
    std::string validationFolder_path = "/Users/sonny/Documents/S&P100/8stocks/targetTime_" + targetTime + "/SAA_kNN";
    std::vector<double> x_SAA = oneStagePiecewiseLP_SAA(folder_path);
    std::cout << "Estimate solution by SAA is " << std::endl;
    for (int x_index = 0; x_index < x_SAA.size(); ++x_index) {
        std::cout << "x[" << x_index << "] = " << x_SAA[x_index] << std::endl;
    }
    // get the objective cost
    double trueCost = 0;
    trueCost = x_SAA[0];
    int base_index = 1;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        trueCost -= exchangeRate_lambda * W_trueResponse[investment_index] * x_SAA[base_index +investment_index];
    }
    // obtain the value in the piecewise linear function
    double PL_value = 0;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        PL_value += -W_trueResponse[investment_index] * x_SAA[base_index + investment_index];
    }
    PL_value += initialValue - x_SAA[0];
    if (PL_value > 0) {
        trueCost += (1.0 / (1.0 - cVaRRisk_alpha)) * PL_value;
    }
    std::cout << "True Cost (by using the true outcome): " << trueCost << std::endl;
    // get the lowest cost of target day
    double lowestTrueCost = oneStagePortfolioOptimizationValidation(folder_path, x_SAA, W_trueResponse, cVaRRisk_alpha, exchangeRate_lambda, initialValue, num_investments);
    validationResult result = oneStagePiecewiseLP_SAA_validation(validationFolder_path, x_SAA);
    std::time_t currTime = std::time(nullptr);
    const char* outputResults_path_const = outputResults_path.c_str();
    std::fstream writeFile;
    writeFile.open(outputResults_path_const,std::fstream::app);
    writeFile << "***************************************************\n";
    writeFile << "Sample Average Approximation has been used\n";
    writeFile << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution: ";
    long x_size = x_SAA.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_SAA[index] << ", ";
    }
    writeFile << x_SAA[x_size - 1] << "\n";
    writeFile << "True Cost (by using the true outcome): " << trueCost << std::endl;
    std::cout << "True Cost (by using the true outcome): " << trueCost << std::endl;
    writeFile << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    std::cout << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "Validation Set Generated by Getting the k Nearest Neighbors of the whole Training Set" << std::endl;
    writeFile << "Mean of cost: " << result.mean << std::endl;
    writeFile << "Variance of cost: " << result.variance << std::endl;
    writeFile << result.alpha << "% CI Lower Bound: " << result.CI_lower << std::endl;
    writeFile << result.alpha << "% CI Upper Bound: " << result.CI_upper << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "***************************************************\n";
    writeFile.close();
}

void oneStagePortfolioOtimization_SAAV2(std::string targetTime, double W_trueResponse[], int num_investments) {
    double initialValue = 100000;
    double cVaRRisk_alpha = 0.1;
    double exchangeRate_lambda = 0.2;
    std::string folder_path = "/Users/sonny/Documents/S&P100/8stocks2.0/targetTime_" + targetTime + "/SAA";
    std::string outputResults_path = folder_path + "/validationResults.txt";
    std::string validationFolder_path = "/Users/sonny/Documents/S&P100/8stocks2.0/targetTime_" + targetTime + "/SAA_kNN";
    std::vector<double> x_SAA = oneStagePiecewiseLP_SAA(folder_path);
    std::cout << "Estimate solution by SAA is " << std::endl;
    for (int x_index = 0; x_index < x_SAA.size(); ++x_index) {
        std::cout << "x[" << x_index << "] = " << x_SAA[x_index] << std::endl;
    }
    // get the objective cost
    double trueCost = 0;
    trueCost = x_SAA[0];
    int base_index = 1;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        trueCost -= exchangeRate_lambda * W_trueResponse[investment_index] * x_SAA[base_index +investment_index];
    }
    // obtain the value in the piecewise linear function
    double PL_value = 0;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        PL_value += -W_trueResponse[investment_index] * x_SAA[base_index + investment_index];
    }
    PL_value += initialValue - x_SAA[0];
    if (PL_value > 0) {
        trueCost += (1.0 / (1.0 - cVaRRisk_alpha)) * PL_value;
    }
    std::cout << "True Cost (by using the true outcome): " << trueCost << std::endl;
    // get the lowest cost of target day
    double lowestTrueCost = oneStagePortfolioOptimizationValidation(folder_path, x_SAA, W_trueResponse, cVaRRisk_alpha, exchangeRate_lambda, initialValue, num_investments);
    validationResult result = oneStagePiecewiseLP_SAA_validation(validationFolder_path, x_SAA);
    std::time_t currTime = std::time(nullptr);
    const char* outputResults_path_const = outputResults_path.c_str();
    std::fstream writeFile;
    writeFile.open(outputResults_path_const,std::fstream::app);
    writeFile << "***************************************************\n";
    writeFile << "Sample Average Approximation has been used\n";
    writeFile << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution: ";
    long x_size = x_SAA.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_SAA[index] << ", ";
    }
    writeFile << x_SAA[x_size - 1] << "\n";
    writeFile << "True Cost (by using the true outcome): " << trueCost << std::endl;
    std::cout << "True Cost (by using the true outcome): " << trueCost << std::endl;
    writeFile << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    std::cout << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "Validation Set Generated by Getting the k Nearest Neighbors of the whole Training Set" << std::endl;
    writeFile << "Mean of cost: " << result.mean << std::endl;
    writeFile << "Variance of cost: " << result.variance << std::endl;
    writeFile << result.alpha << "% CI Lower Bound: " << result.CI_lower << std::endl;
    writeFile << result.alpha << "% CI Upper Bound: " << result.CI_upper << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "***************************************************\n";
    writeFile.close();
}

void robustOneStagePortfolioOtimization_SSubgradient(std::vector<double> x, int maxIterations) {
    double M = 307.6386;
    double lowerBound = 0;
    double upperBound = 200;
    double Dx = Dx_estimateRectangle(x, 0, 200);
    //std::string caseNumber = "2";
    std::string targetTime = "2015";
    // true response
    double W_trueResponse[] = {1, 139.75705, 47.207386, 119.690994, 141.210007, 235.925293, 59.700073, 85.5019, 164.313843};
    std::cout << "Dx: " << Dx << std::endl;
    int m = 5;
    //std::string folder_path = "/Users/sonny/Documents/S&P100/case" + caseNumber + "/" + method;
    std::string folder_path = "/Users/sonny/Documents/S&P100/8stocks/targetTime_" + targetTime + "/stochasticSubgradient";
    std::string validationFolder_path = "/Users/sonny/Documents/S&P100/8stocks/targetTime_" + targetTime + "/SAA_kNN";
    // solve by using robust NSQG
    std::vector<double> x_candidate = robustOneStagePiecewiseLP_outputResults(folder_path, maxIterations, x, m, M, lowerBound, upperBound);
    int num_investments = (sizeof(W_trueResponse) / sizeof(W_trueResponse[0]));
    std::cout << num_investments << std::endl;
    double initialValue = 100000;
    double cVaRRisk_alpha = 0.1;
    double exchangeRate_lambda = 0.2;
    std::cout << "Start Measuring the Quality of the Solution" << std::endl;
    // get the lowest cost of target day
    double lowestTrueCost = oneStagePortfolioOptimizationValidation(folder_path, x_candidate, W_trueResponse, cVaRRisk_alpha, exchangeRate_lambda, initialValue, num_investments);
    // get the objective cost
    double trueCost = 0;
    trueCost = x_candidate[0];
    int base_index = 1;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        trueCost -= exchangeRate_lambda * W_trueResponse[investment_index] * x_candidate[base_index +investment_index];
    }
    // obtain the value in the piecewise linear function
    double PL_value = 0;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        PL_value += -W_trueResponse[investment_index] * x_candidate[base_index + investment_index];
    }
    PL_value += initialValue - x_candidate[0];
    if (PL_value > 0) {
        trueCost += (1.0 / (1.0 - cVaRRisk_alpha)) * PL_value;
    }
    validationResult result = oneStagePiecewiseLP_SAA_validation(validationFolder_path, x_candidate);
    // write the results of quality of candidate solution
    std::string outputResults_path = folder_path + "/validationResults.txt";
    const char* outputResults_path_const = outputResults_path.c_str();
    std::fstream writeFile;
    writeFile.open(outputResults_path_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "***************************************************\n";
    writeFile << "Robust Stochastic Subgradient Method has been used" << std::endl;
    writeFile << "Number of outer loops: " << maxIterations << std::endl;
    writeFile << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution: ";
    long x_size = x_candidate.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_candidate[index] << ", ";
    }
    writeFile << x_candidate[x_size - 1] << "\n";
    writeFile << "True Cost (by using the true outcome): " << trueCost << std::endl;
    std::cout << "True Cost (by using the true outcome): " << trueCost << std::endl;
    writeFile << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    std::cout << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "Validation Set Generated by Getting the k Nearest Neighbors of the whole Training Set" << std::endl;
    writeFile << "Mean of cost: " << result.mean << std::endl;
    writeFile << "Variance of cost: " << result.variance << std::endl;
    writeFile << result.alpha << "% CI Lower Bound: " << result.CI_lower << std::endl;
    writeFile << result.alpha << "% CI Upper Bound: " << result.CI_upper << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "***************************************************\n";
    writeFile.close();
}

void robustOneStagePortfolioOtimization_SSubgradientV2(std::vector<double> x, int maxIterations, std::string targetTime, double W_trueResponse[], int num_investments) {
    double M = 307.6386;
    double lowerBound = 0;
    double upperBound = 200;
    double Dx = Dx_estimateRectangle(x, 0, 200);
    std::cout << "Dx: " << Dx << std::endl;
    int m = 5;
    //std::string folder_path = "/Users/sonny/Documents/S&P100/case" + caseNumber + "/" + method;
    std::string folder_path = "/Users/sonny/Documents/S&P100/8stocks2.0/targetTime_" + targetTime + "/stochasticSubgradient";
    std::string validationFolder_path = "/Users/sonny/Documents/S&P100/8stocks2.0/targetTime_" + targetTime + "/SAA_kNN";
    // solve by using robust NSQG
    std::vector<double> x_candidate = robustOneStagePiecewiseLP_outputResults(folder_path, maxIterations, x, m, M, lowerBound, upperBound);
    std::cout << num_investments << std::endl;
    double initialValue = 100000;
    double cVaRRisk_alpha = 0.1;
    double exchangeRate_lambda = 0.2;
    std::cout << "Start Measuring the Quality of the Solution" << std::endl;
    // get the lowest cost of target day
    double lowestTrueCost = oneStagePortfolioOptimizationValidation(folder_path, x_candidate, W_trueResponse, cVaRRisk_alpha, exchangeRate_lambda, initialValue, num_investments);
    // get the objective cost
    double trueCost = 0;
    trueCost = x_candidate[0];
    int base_index = 1;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        trueCost -= exchangeRate_lambda * W_trueResponse[investment_index] * x_candidate[base_index +investment_index];
    }
    // obtain the value in the piecewise linear function
    double PL_value = 0;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        PL_value += -W_trueResponse[investment_index] * x_candidate[base_index + investment_index];
    }
    PL_value += initialValue - x_candidate[0];
    if (PL_value > 0) {
        trueCost += (1.0 / (1.0 - cVaRRisk_alpha)) * PL_value;
    }
    validationResult result = oneStagePiecewiseLP_SAA_validation(validationFolder_path, x_candidate);
    // write the results of quality of candidate solution
    std::string outputResults_path = folder_path + "/validationResults.txt";
    const char* outputResults_path_const = outputResults_path.c_str();
    std::fstream writeFile;
    writeFile.open(outputResults_path_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "***************************************************\n";
    writeFile << "Robust Stochastic Subgradient Method has been used" << std::endl;
    writeFile << "Number of outer loops: " << maxIterations << std::endl;
    writeFile << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution: ";
    long x_size = x_candidate.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_candidate[index] << ", ";
    }
    writeFile << x_candidate[x_size - 1] << "\n";
    writeFile << "True Cost (by using the true outcome): " << trueCost << std::endl;
    std::cout << "True Cost (by using the true outcome): " << trueCost << std::endl;
    writeFile << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    std::cout << "Lowest true cost (by solving the portfolio optimization problem with true outcome): " << lowestTrueCost << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "Validation Set Generated by Getting the k Nearest Neighbors of the whole Training Set" << std::endl;
    writeFile << "Mean of cost: " << result.mean << std::endl;
    writeFile << "Variance of cost: " << result.variance << std::endl;
    writeFile << result.alpha << "% CI Lower Bound: " << result.CI_lower << std::endl;
    writeFile << result.alpha << "% CI Upper Bound: " << result.CI_upper << std::endl;
    writeFile << "--------------------------------------------------\n";
    writeFile << "***************************************************\n";
    writeFile.close();
}

void NSQG_user_interface(){
    // parameters
    int methodNumber = 0;
    char ifRobust = 'N';
    double initStepsize = 0.5;
    std::string initialPointInput;
    std::string folderPath;
    std::cout << "Please enter the folder path for database: " << std::endl;
    std::cin >> folderPath;
    std::cout << "Current available methods are listed below: " << std::endl;
    std::cout << "(1) kNN\n";
    std::cout << "(2) naive kernel\n";
    std::cout << "(3) quartic kernel\n";
    std::cout << "(4) Epanechnikov kernel\n";
    std::cout << "Please specify the method by entering the correponding number" << std::endl;
    std::cout << "Method number (e.g. 1): ";
    std::cin >> methodNumber;
    std::cout << "Is robust version used? (Y/N): ";
    std::cin >> ifRobust;
    std::cout << "Please enter the initial stepsize: ";
    std::cin >> initStepsize;
    std::cout << "Please enter the initial point: ";
    std::cin >> initialPointInput;
    // extract the initial point
    std::stringstream ss1(initialPointInput);
    std::string line1;
    std::vector<double> initialPoint;
    std::cout << "The initial point that you enter is: ";
    while (getline(ss1, line1, ',')) {
        double entry;
        std::stringstream ss2(line1);
        ss2 >> entry;
        initialPoint.push_back(entry);
        std::cout << entry << " ";
    }
    std::cout << std::endl;
    if (ifRobust == 'N') {
        std::string trainDB_path;
        switch (methodNumber) {
            case 1:
                std::cout << "SQG-kNN is chosen.\n";
                trainDB_path = folderPath + "/kNN";
                break;
            case 2:
                std::cout << "SQG-naiveKernel is chosen.\n";
                trainDB_path = folderPath + "/naiveKernel";
                break;
            case 3:
                std::cout << "SQG-quarticKernel is chosen.\n";
                trainDB_path = folderPath + "/quarticKernel";
                break;
            case 4:
                std::cout << "SQG-EpanechnikovKernel is chosen.\n";
                trainDB_path = folderPath + "/EpanechnikovKernel";
                break;
            default:
                break;
        }
        int maxIterations = 5; // 5 15 30 50 75
        std::cout << "Please enter the maximum number of iterations: " << std::endl;
        std::cin >> maxIterations;
        std::cout << "######################################################\n";
        std::cout << "Initiate NSQG\n";
        std::vector<double> x_candidate = twoStageLP_random_b_outputResults(trainDB_path, maxIterations, initStepsize, initialPoint);
        // validation
        std::string validationFolder_path = folderPath + "/trueDist";
        std::string empiricalValidationFolder_path = folderPath + "/empiricalDist";
        // estimate the quality of solution based on the true validation set
        validationResult trueResult = twoStageLP_validation_outputResultsV2(validationFolder_path, x_candidate);
        // estimate the quality of solution based on the kNN validation set
        validationResult empiricalValidationResult = twoStageLP_validation_outputResultsV2(empiricalValidationFolder_path, x_candidate);
        // write the results of quality of candidate solution
        std::string outputResults_path = folderPath + "/validationResults.txt";
        const char* outputResults_path_const = outputResults_path.c_str();
        std::fstream writeFile;
        writeFile.open(outputResults_path_const,std::fstream::app);
        // current time
        std::time_t currTime = std::time(nullptr);
        std::cout << "**********************************************************************\n";
        writeFile << "**********************************************************************\n";
        switch (methodNumber) {
            case 1:
                writeFile << "Nonparametric SQG-kNN has been used" << std::endl;
                break;
            case 2:
                writeFile << "Nonparametric SQG-naiveKernel has been used" << std::endl;
                break;
            case 3:
                writeFile << "Nonparametric SQG-quarticKernel has been used" << std::endl;
                break;
            case 4:
                writeFile << "Nonparametric SQG-EpanechnikovKernel has been used" << std::endl;
                break;
            default:
                break;
        }
        writeFile << "Initial Point        : ";
        std::cout << "Initial Point        : ";
        long x_size = x_candidate.size();
        for (int index = 0; index < x_size - 1; ++index) {
            writeFile << initialPoint[index] << ", ";
            std::cout << initialPoint[index] << ", ";
        }
        writeFile << initialPoint[x_size - 1] << "\n";
        std::cout << initialPoint[x_size - 1] << "\n";
        writeFile << "Initial Step Size    : " << initStepsize << "\n";
        std::cout << "Initial Step Size    : " << initStepsize << "\n";
        writeFile << "Number of Iterations : " << maxIterations << "\n";
        std::cout << "Number of Iterations : " << maxIterations << "\n";
        writeFile << "Estimating the quality of candidate solution\n";
        std::cout << "Estimating the quality of candidate solution\n";
        writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
        std::cout << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
        writeFile << "Candidate solution   : ";
        std::cout << "Candidate solution   : ";
        for (int index = 0; index < x_size - 1; ++index) {
            writeFile << x_candidate[index] << ", ";
            std::cout << x_candidate[index] << ", ";
        }
        writeFile << x_candidate[x_size - 1] << "\n";
        std::cout << x_candidate[x_size - 1] << "\n";
        writeFile << "========================True Validation Set========================\n";
        std::cout << "========================True Validation Set========================\n";
        writeFile << "Number of data points: " << trueResult.num_dataPoint << std::endl;
        std::cout << "Number of data points: " << trueResult.num_dataPoint << std::endl;
        writeFile << "Average cost         : " << trueResult.mean << std::endl;
        std::cout << "Average cost         : " << trueResult.mean << std::endl;
        writeFile << "Variance             : " << trueResult.variance << std::endl;
        std::cout << "Variance             : " << trueResult.variance << std::endl;
        writeFile << "Error in estimating the expected cost: " << sqrt(trueResult.variance / ((double)trueResult.num_dataPoint)) << std::endl;
        std::cout << "Error in estimating the expected cost: " << sqrt(trueResult.variance / ((double)trueResult.num_dataPoint)) << std::endl;
        writeFile << trueResult.alpha << "% CI Lower Bound: " << trueResult.CI_lower << std::endl;
        std::cout << trueResult.alpha << "% CI Lower Bound: " << trueResult.CI_lower << std::endl;
        writeFile << trueResult.alpha << "% CI Upper Bound: " << trueResult.CI_upper << std::endl;
        std::cout << trueResult.alpha << "% CI Upper Bound: " << trueResult.CI_upper << std::endl;
        writeFile << "========================kNN Validation Set=========================\n";
        std::cout << "========================kNN Validation Set=========================\n";
        writeFile << "Number of data points: " << empiricalValidationResult.num_dataPoint << std::endl;
        std::cout << "Number of data points: " << empiricalValidationResult.num_dataPoint << std::endl;
        writeFile << "Average cost         : " << empiricalValidationResult.mean << std::endl;
        std::cout << "Average cost         : " << empiricalValidationResult.mean << std::endl;
        writeFile << "Variance             : " << empiricalValidationResult.variance << std::endl;
        std::cout << "Variance             : " << empiricalValidationResult.variance << std::endl;
        writeFile << "Error in estimating the expected cost: " << sqrt(empiricalValidationResult.variance / ((double)empiricalValidationResult.num_dataPoint)) << std::endl;
        std::cout << "Error in estimating the expected cost: " << sqrt(empiricalValidationResult.variance / ((double)empiricalValidationResult.num_dataPoint)) << std::endl;
        writeFile << empiricalValidationResult.alpha << "% CI Lower Bound: " << empiricalValidationResult.CI_lower << std::endl;
        std::cout << empiricalValidationResult.alpha << "% CI Lower Bound: " << empiricalValidationResult.CI_lower << std::endl;
        writeFile << trueResult.alpha << "% CI Upper Bound: " << empiricalValidationResult.CI_upper << std::endl;
        std::cout << empiricalValidationResult.alpha << "% CI Upper Bound: " << empiricalValidationResult.CI_upper << std::endl;
        writeFile << "**********************************************************************\n";
        std::cout << "**********************************************************************\n";
        writeFile.close();
    }
    else if (ifRobust == 'Y') {
        std::string trainDB_path;
        switch (methodNumber) {
            case 1:
                std::cout << "Robust SQG-kNN is chosen.\n";
                trainDB_path = folderPath + "/kNN";
                break;
            case 2:
                std::cout << "Robust SQG-naiveKernel is chosen.\n";
                trainDB_path = folderPath + "/naiveKernel";
                break;
            case 3:
                std::cout << "Robust SQG-quarticKernel is chosen.\n";
                trainDB_path = folderPath + "/quarticKernel";
                break;
            case 4:
                std::cout << "Robust SQG-EpanechnikovKernel is chosen.\n";
                trainDB_path = folderPath + "/EpanechnikovKernel";
            default:
                break;
        }
        // parameters for the robust NSQG
        int maxOuterLoops = 5; // 1 2 3 4 5
        int m = 5;
        double M = 30;
        std::cout << "Please enter the max number of outer loops: ";
        std::cin >> maxOuterLoops;
        std::cout << "Please enter the increment in the inner loops (e.g. enter 5): ";
        std::cin >> m;
        std::cout << "Please enter the upper bound of the subgradient (e.g. enter 30): ";
        std::cin >> M;
        std::cout << "######################################################\n";
        std::cout << "Initiate Robust NSQG\n";
        double Dx = Dx_estimateRectangle(initialPoint, 0, 20);
        std::vector<double> x_candidate = robustTwoStageLP_random_b_outputResults(trainDB_path, maxOuterLoops, initialPoint, Dx, m, M);
        // validation
        std::string validationFolder_path = folderPath + "/trueDist";
        std::string empiricalValidationFolder_path = folderPath + "/empiricalDist";
        // estimate the quality of solution based on the true validation set
        validationResult trueResult = twoStageLP_validation_outputResultsV2(validationFolder_path, x_candidate);
        // estimate the quality of solution based on the kNN validation set
        validationResult empiricalValidationResult = twoStageLP_validation_outputResultsV2(empiricalValidationFolder_path, x_candidate);
        // write the results of quality of candidate solution
        std::string outputResults_path = folderPath + "/validationResults.txt";
        const char* outputResults_path_const = outputResults_path.c_str();
        std::fstream writeFile;
        writeFile.open(outputResults_path_const,std::fstream::app);
        // current time
        std::time_t currTime = std::time(nullptr);
        std::cout << "**********************************************************************\n";
        writeFile << "**********************************************************************\n";
        switch (methodNumber) {
            case 1:
                writeFile << "Robust Nonparametric SQG-kNN has been used" << std::endl;
                break;
            case 2:
                writeFile << "Robust Nonparametric SQG-naiveKernel has been used" << std::endl;
                break;
            case 3:
                writeFile << "Robust Nonparametric SQG-quarticKernel has been used" << std::endl;
                break;
            case 4:
                writeFile << "Robust Nonparametric SQG-EpanechnikovKernel has been used" << std::endl;
                break;
            default:
                break;
        }
        writeFile << "Initial Point        : ";
        std::cout << "Initial Point        : ";
        long x_size = x_candidate.size();
        for (int index = 0; index < x_size - 1; ++index) {
            writeFile << initialPoint[index] << ", ";
            std::cout << initialPoint[index] << ", ";
        }
        writeFile << initialPoint[x_size - 1] << "\n";
        std::cout << initialPoint[x_size - 1] << "\n";
        writeFile << "Number of Iterations (Outer loop): " << maxOuterLoops << "\n";
        std::cout << "Number of Iterations (Outer loop): " << maxOuterLoops << "\n";
        writeFile << "Dx                   : " << Dx << "\n";
        std::cout << "Dx                   : " << Dx << "\n";
        writeFile << "M                    : " << M << "\n";
        std::cout << "M                    : " << M << "\n";
        writeFile << "m                    : " << m << "\n";
        std::cout << "m                    : " << m << "\n";
        writeFile << "Estimating the quality of candidate solution\n";
        std::cout << "Estimating the quality of candidate solution\n";
        writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
        std::cout << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
        writeFile << "Candidate solution   : ";
        std::cout << "Candidate solution   : ";
        for (int index = 0; index < x_size - 1; ++index) {
            writeFile << x_candidate[index] << ", ";
            std::cout << x_candidate[index] << ", ";
        }
        writeFile << x_candidate[x_size - 1] << "\n";
        std::cout << x_candidate[x_size - 1] << "\n";
        writeFile << "========================True Validation Set========================\n";
        std::cout << "========================True Validation Set========================\n";
        writeFile << "Number of data points: " << trueResult.num_dataPoint << std::endl;
        std::cout << "Number of data points: " << trueResult.num_dataPoint << std::endl;
        writeFile << "Average cost         : " << trueResult.mean << std::endl;
        std::cout << "Average cost         : " << trueResult.mean << std::endl;
        writeFile << "Variance             : " << trueResult.variance << std::endl;
        std::cout << "Variance             : " << trueResult.variance << std::endl;
        writeFile << "Error in estimating the expected cost: " << sqrt(trueResult.variance / ((double)trueResult.num_dataPoint)) << std::endl;
        std::cout << "Error in estimating the expected cost: " << sqrt(trueResult.variance / ((double)trueResult.num_dataPoint)) << std::endl;
        writeFile << trueResult.alpha << "% CI Lower Bound: " << trueResult.CI_lower << std::endl;
        std::cout << trueResult.alpha << "% CI Lower Bound: " << trueResult.CI_lower << std::endl;
        writeFile << trueResult.alpha << "% CI Upper Bound: " << trueResult.CI_upper << std::endl;
        std::cout << trueResult.alpha << "% CI Upper Bound: " << trueResult.CI_upper << std::endl;
        writeFile << "========================kNN Validation Set=========================\n";
        std::cout << "========================kNN Validation Set=========================\n";
        writeFile << "Number of data points: " << empiricalValidationResult.num_dataPoint << std::endl;
        std::cout << "Number of data points: " << empiricalValidationResult.num_dataPoint << std::endl;
        writeFile << "Average cost         : " << empiricalValidationResult.mean << std::endl;
        std::cout << "Average cost         : " << empiricalValidationResult.mean << std::endl;
        writeFile << "Variance             : " << empiricalValidationResult.variance << std::endl;
        std::cout << "Variance             : " << empiricalValidationResult.variance << std::endl;
        writeFile << "Error in estimating the expected cost: " << sqrt(empiricalValidationResult.variance / ((double)empiricalValidationResult.num_dataPoint)) << std::endl;
        std::cout << "Error in estimating the expected cost: " << sqrt(empiricalValidationResult.variance / ((double)empiricalValidationResult.num_dataPoint)) << std::endl;
        writeFile << empiricalValidationResult.alpha << "% CI Lower Bound: " << empiricalValidationResult.CI_lower << std::endl;
        std::cout << empiricalValidationResult.alpha << "% CI Lower Bound: " << empiricalValidationResult.CI_lower << std::endl;
        writeFile << trueResult.alpha << "% CI Upper Bound: " << empiricalValidationResult.CI_upper << std::endl;
        std::cout << empiricalValidationResult.alpha << "% CI Upper Bound: " << empiricalValidationResult.CI_upper << std::endl;
        writeFile << "**********************************************************************\n";
        std::cout << "**********************************************************************\n";
        writeFile.close();
    }
}

int main(int argc, const char * argv[]) {
    NSQG_user_interface();
    /*
    //inputOneStageModelTest();
    //testStochasticOneStagePiecewiseLP();
    // test on NSQG
    std::vector<double> x(4,4.0);
    int caseNumber = 1;
    int validationCaseNumber = 1 ;
    int kNNValidationCaseNumber = 1;
    int experimentNumber = 4; // 2, 3, 4
    //std::string method = "kNN";
    //std::string method = "naiveKernel";
    //std::string method = "quarticKernel";
    //std::string method = "EpanechnikovKernel";
    //std::string method = "kNN3";
    std::string method = "naiveKernel3";
    //std::string method = "quarticKernel3";
    //std::string method = "EpanechnikovKernel3";
    int maxIterations = 5; // 5 15 30 50 75
    int maxOuterLoops = 5; // 1 2 3 4 5
    double initStepsize = 0.8;
    //NSQG_tssV2_2(caseNumber, method, x, maxIterations, initStepsize, validationCaseNumber, kNNValidationCaseNumber, experimentNumber);
    //SS_tssV2_2(caseNumber, x, maxIterations, initStepsize, validationCaseNumber, kNNValidationCaseNumber, experimentNumber);
    //robustNSQG_tssV2_2(caseNumber, method, x, maxOuterLoops, validationCaseNumber, kNNValidationCaseNumber,experimentNumber);
    //robustSS_tssV2_2(caseNumber, x, maxOuterLoops, validationCaseNumber, kNNValidationCaseNumber,experimentNumber);
    //NSQG_tssV2(caseNumber,method,x,maxIterations,initStepsize);
    //SAA_tsskNNValidation(kNNValidationCaseNumber);
    SAA_tssTrueValidation(validationCaseNumber);
     */
    //===============================================================================
    // Estimate error
    //double x0[] = {7.46378, 10.126, 11.6792, 9.74133};
    /*
    double x0[] = {4, 4, 4, 4};
    std::vector<double> x_candidtate(4,0.0);
    for (int x_index = 0; x_index < 4; ++x_index) {
        x_candidtate[x_index] = x0[x_index];
    }
    int caseNumber = 28; // 1 2 ... 30
    int experimentNumber = 4; // 2, 3, 4
    int dataset_index = 90; //0 to 99
    int num_cases = 60;
    //std::string method = "kNN";
    //std::string method = "naiveKernel";
    //std::string method = "quarticKernel";
    //std::string method = "EpanechnikovKernel";
    //std::string method = "kNN3";
    std::string method = "naiveKernel3";
    //std::string method = "quarticKernel3";
    //std::string method = "EpanechnikovKernel3";
    std::cout << "Error Calculation\n";
    //tss_errorCalculation(caseNumber, method, x_candidtate, experimentNumber, dataset_index);
    //tss_errorCalculationV2(experimentNumber, method, x_candidtate, dataset_index);
    tss_errorCalculationV3(experimentNumber, method, x_candidtate, dataset_index, num_cases);
     */
    //===============================================================================
    // test on robust NSQG
    /*
     std::vector<double> x(4,4);
     int caseNumber = 3;
     //std::string method = "kNN";
     std::string method = "naiveKernel";
     //std::string method = "quarticKernel";
     //std::string method = "EpanechnikovKernel";
     int maxIterations = 1;
     //robustNSQG_tss(caseNumber, method, x, maxIterations);
     robustNSQG_tssV2(caseNumber, method, x, maxIterations);
     //projectionTest();
     */
    // test on NSQG for portfolio Optimizaion
    /*
    std::vector<double> x(28,0);
    x[1] = 100000;
    double initialStepsize = 10;
    int maxIterations = 800;
    int maxIterations_robust = 17;
    //std::string caseNumber = "2";
    std::string targetTime = "2022";
    // true response
    double W_trueResponse[] = {1, 144.558197, 47.381447, 122.247871, 146.160004, 237.43187, 60.964828, 85.732979, 167.02002};
    int num_investment = (sizeof(W_trueResponse) / sizeof(W_trueResponse[0]));
    //std::string method = "kNN";
    //std::string method = "naiveKernel";
    //std::string method = "EpanechnikovKernel";
    std::string method = "quarticKernel";
    //NSQG_portfolioOptimization(method, x, maxIterations, initialStepsize);
    //NSQG_portfolioOptimizationV2(method, x, maxIterations, initialStepsize, targetTime, W_trueResponse, num_investment);
    //robustNSQG_portfolioOptimization(method, x, maxIterations_robust);
    //robustNSQG_portfolioOptimizationV2(method, x, maxIterations_robust, targetTime, W_trueResponse, num_investment);
    //robustNSQG_portfolioOptimizationDEBUG(method, x, maxIterations_robust);
    //oneStageValidationTest();
    //oneStagePortfolioOtimization_SAA();
    oneStagePortfolioOtimization_SAAV2(targetTime, W_trueResponse, num_investment);
    //robustOneStagePortfolioOtimization_SSubgradient(x, maxIterations_robust);
    robustOneStagePortfolioOtimization_SSubgradientV2(x, maxIterations_robust, targetTime, W_trueResponse, num_investment);
     */
    return 0;
}
