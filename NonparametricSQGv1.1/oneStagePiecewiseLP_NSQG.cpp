//
//  oneStagePiecewiseLP_NSQG.cpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/24/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#include "oneStagePiecewiseLP_NSQG.hpp"

// One Stage Piecewise Linear Programming
std::vector<double> oneStagePiecewiseLP_outputResults(const std::string& folder_path, int max_iterations, double initial_stepsize, std::vector<double> x_init) {
    // STEP 1: INITIALIZATION
    // create directory paths for database and model
    std::string DB_path = folder_path + "/DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults.txt";
    // create database
    std::vector<std::vector<dataPoint>> rawObjectiveDB = readNonparametricDB(DB_path);
    // create model structure
    oneStageParameters parameters = readOneStageParameters(model_path);
    // size of x
    long x_size = parameters.x.size();
    // create sto object
    oneStagePiecewiseLinearObjective parameters_objective = readStochasticOneStagePiecewiseLP(sto_path);
    // merge database and parameters of objeective
    oneStagePiecewiseLinearObjectiveDB objectiveDB = mergeDB_randomObjective(rawObjectiveDB, parameters_objective, x_size);
    // STEP 2: SOLVE
    std::vector<double> x_est = oneStagePiecewiseLP_solver(parameters, objectiveDB, initial_stepsize, max_iterations, x_init, resultsOutput_path);
    return x_est; // return the estimated solution
} // end oneStagePiecewiseLP_outputResults


// One Stage Piecewise Linear Programming Solved by Robust Nonparametric SQG
std::vector<double> robustOneStagePiecewiseLP_outputResults(const std::string& folder_path, int max_iterations, std::vector<double> x_init, int m, double M, double lowerBound, double upperBound) {
    // STEP 1: INITIALIZATION
    // create directory paths for database and model
    std::string DB_path = folder_path + "/DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(Robust).txt";
    // create database
    std::vector<std::vector<dataPoint>> rawObjectiveDB = readNonparametricDB(DB_path);
    // create model structure
    oneStageParameters parameters = readOneStageParameters(model_path);
    // size of x
    long x_size = parameters.x.size();
    // create sto object
    oneStagePiecewiseLinearObjective parameters_objective = readStochasticOneStagePiecewiseLP(sto_path);
    // calculate Dx
    double Dx = DxEstimate_finer(x_init, parameters, lowerBound, upperBound);
    std::cout << "Dx(finer): " << Dx << std::endl;
    // merge database and parameters of objeective
    oneStagePiecewiseLinearObjectiveDB objectiveDB = mergeDB_randomObjective(rawObjectiveDB, parameters_objective, x_size);
    // STEP 2: SOLVE
    std::vector<double> x_est = robustOneStagePiecewiseLP_solver(parameters, objectiveDB, max_iterations, x_init, Dx, m, M, resultsOutput_path);
    return x_est; // return the estimated solution
} // end robustOneStagePiecewiseLP_outputResults

// for debugging
std::vector<double> robustOneStagePiecewiseLP_outputResultsDEBUG(const std::string& folder_path, int max_iterations, std::vector<double> x_init, int m, double M, double lowerBound, double upperBound) {
    // STEP 1: INITIALIZATION
    // create directory paths for database and model
    std::string DB_path = folder_path + "/DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string resultsOutput_path = folder_path + "/DEBUGcomputationalResults(Robust).txt";
    // create database
    std::vector<std::vector<dataPoint>> rawObjectiveDB = readNonparametricDB(DB_path);
    // create model structure
    oneStageParameters parameters = readOneStageParameters(model_path);
    // size of x
    long x_size = parameters.x.size();
    // create sto object
    oneStagePiecewiseLinearObjective parameters_objective = readStochasticOneStagePiecewiseLP(sto_path);
    // calculate Dx
    double Dx = DxEstimate_finer(x_init, parameters, lowerBound, upperBound);
    std::cout << "Dx(finer): " << Dx << std::endl;
    // merge database and parameters of objeective
    oneStagePiecewiseLinearObjectiveDB objectiveDB = mergeDB_randomObjective(rawObjectiveDB, parameters_objective, x_size);
    // STEP 2: SOLVE
    std::vector<double> x_est = robustOneStagePiecewiseLP_solverDEBUG(parameters, objectiveDB, max_iterations, x_init, Dx, m, M, resultsOutput_path);
    return x_est; // return the estimated solution
}

// Portforlio Optimization Validation by using true outcome (no prediction)
double oneStagePortfolioOptimizationValidation(const std::string& folder_path, const std::vector<double>& x_est,double W_trueResponse[], double cVaRRisk_alpha, double exchangeRate_lambda, double initialValue, int num_investments) {
    std::cout << "Obtaining Optimal Cost by Knowing the True Outcome" << std::endl;
    double trueCost = 0;
    // create directory path for model
    std::string model_path = folder_path + "/model.txt";
    // size of x
    long x_size = x_est.size();
    std::cout << "number of investments: "<< num_investments << std::endl;
    // read model structure
    oneStageParameters parameters = readOneStageParameters(model_path);
    // number of rows in A
    long A_rowsize = parameters.A.size();
    // create cplex solver environment for getting the lowest possible risk
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray x(env, x_size, -IloInfinity, IloInfinity, ILOFLOAT);
    IloNumVar y(env, 0, IloInfinity, ILOFLOAT);
    mod.add(x);
    mod.add(y);
    // objective function
    IloExpr expr_obj(env);
    expr_obj = x[0];
    int base_index = 1;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        expr_obj += -exchangeRate_lambda * W_trueResponse[investment_index] * x[base_index + investment_index];
    }
    expr_obj += (1.0 / (1.0 - cVaRRisk_alpha)) * y;
    IloObjective obj = IloMinimize(env, expr_obj);
    mod.add(obj);
    // constraints
    // constraint for the piecewise linear function
    IloExpr expr_constraintPL(env);
    expr_constraintPL = y;
    for (int investment_index = 0; investment_index < num_investments; ++investment_index) {
        expr_constraintPL += W_trueResponse[investment_index] * x[base_index + investment_index];
    }
    expr_constraintPL += x[0] - initialValue;
    mod.add(expr_constraintPL >= 0);
    // other constraints
    for (int row_index = 0; row_index < A_rowsize; ++row_index) {
        IloExpr expr_constrant(env);
        for (auto a_iterator : parameters.A[row_index]) {
            int x_index = a_iterator.first;
            double coefficient = a_iterator.second;
            //std::cout << "index of x: " << x_index << std::endl;
            expr_constrant += coefficient * x[x_index];
        }
        // b is stored in the form of hash map, the iterator corresponding to the target key needs to be found
        std::unordered_map<int, double>::const_iterator b_iterator = parameters.b.find(row_index);
        if (b_iterator == parameters.b.end()) { // if b[row_index] is 0
            mod.add(expr_constrant <= 0);
        }
        else { // if b[row_index] is not 0
            mod.add(expr_constrant <= b_iterator->second);
        }
    }
    // create cplex environment
    IloCplex cplex(env);
    cplex.extract(mod);
    //cplex.setOut(env.getNullStream()); // stop the output
    //cplex.exportModel("/Users/sonny/Documents/S&P100/8stocks/targetTime_2013/model01.lp");
    cplex.solve();
    // obtain optimal cost
    for (int index = 0; index < x_size; ++index) {
        std::cout << "x[" << index << "]: " << cplex.getValue(x[index]) << std::endl;
    }
    std::cout << "y: " << cplex.getValue(y) << std::endl;
    trueCost = cplex.getObjValue();
    env.end();
    return trueCost;
}

// One Stage Piecewise Linear Programming solved by SAA
std::vector<double> oneStagePiecewiseLP_SAA(const std::string& folder_path) {
    // STEP 1: INITIALIZATION
    // create directory paths for database and model
    std::string DB_path = folder_path + "/DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(SAA).txt";
    // create database
    std::vector<std::vector<dataPoint>> rawObjectiveDB = readNonparametricDB(DB_path);
    // create model structure
    oneStageParameters parameters = readOneStageParameters(model_path);
    // number of scenarios
    long num_scenarios = rawObjectiveDB[0].size();
    double fraction = 1.0 / ((double) num_scenarios);
    // size of x
    long x_size = parameters.x.size();
    // number of rows in A
    long A_rowsize = parameters.A.size();
    // create sto object
    oneStagePiecewiseLinearObjective parameters_objective = readStochasticOneStagePiecewiseLP(sto_path);
    // merge database and parameters of objeective
    oneStagePiecewiseLinearObjectiveDB objectiveDB = mergeDB_randomObjective(rawObjectiveDB, parameters_objective, x_size);
    // creat SAA problem
    IloEnv env;
    IloModel mod(env);
    // decision variables
    IloNumVarArray x(env, x_size, -IloInfinity, IloInfinity, ILOFLOAT);
    IloNumVarArray y(env, num_scenarios, -IloInfinity, IloInfinity, ILOFLOAT);
    mod.add(x);
    mod.add(y);
    // objective function
    IloExpr expr_obj(env);
    for (int scenario_index = 0; scenario_index < num_scenarios; ++scenario_index) {
        // c1^T x
        for (auto x_iterator : objectiveDB.c1[0][scenario_index]) {
            int x_index = x_iterator.first; // index of x
            double coefficient = x_iterator.second; // coeffcient for x[x_index]
            expr_obj += fraction * coefficient * x[x_index];
        }
        // e1
        expr_obj += fraction * objectiveDB.e1[0][scenario_index];
        // piecewise linear function
        expr_obj += fraction * objectiveDB.cPiecewise * y[scenario_index];
    }
    IloObjective obj = IloMinimize(env, expr_obj);
    mod.add(obj);
    // constraints
    // constraints for piecewise linear functions
    for (int scenario_index = 0; scenario_index < num_scenarios; ++scenario_index) {
        IloExpr expr_PLconstraint1(env);
        expr_PLconstraint1 = y[scenario_index];
        // c2^T x
        for (auto x_iterator : objectiveDB.c2[0][scenario_index]) {
            int x_index = x_iterator.first; // index of x
            double coefficient = x_iterator.second; // coefficient for x[x_index]
            expr_PLconstraint1 += -coefficient * x[x_index];
        }
        // e2
        expr_PLconstraint1 -= objectiveDB.e2[0][scenario_index];
        mod.add(expr_PLconstraint1 >= 0);
        IloExpr expr_PLconstraint2(env);
        expr_PLconstraint2 = y[scenario_index];
        // c3^T x
        for (auto x_iterator : objectiveDB.c3[0][scenario_index]) {
            int x_index = x_iterator.first; // index of x
            double coefficient = x_iterator.second; // coefficient for x[x_index]
            expr_PLconstraint2 += -coefficient * x[x_index];
        }
        // e3
        expr_PLconstraint2 -= objectiveDB.e3[0][scenario_index];
        mod.add(expr_PLconstraint2 >= 0);
    }
    // Ax <= b
    for (int row_index = 0; row_index < A_rowsize; ++row_index) {
        IloExpr expr_constrant(env);
        for (auto a_iterator : parameters.A[row_index]) {
            int x_index = a_iterator.first;
            double coefficient = a_iterator.second;
            //std::cout << "index of x: " << x_index << std::endl;
            expr_constrant += coefficient * x[x_index];
        }
        // b is stored in the form of hash map, the iterator corresponding to the target key needs to be found
        std::unordered_map<int, double>::const_iterator b_iterator = parameters.b.find(row_index);
        if (b_iterator == parameters.b.end()) { // if b[row_index] is 0
            mod.add(expr_constrant <= 0);
        }
        else { // if b[row_index] is not 0
            mod.add(expr_constrant <= b_iterator->second);
        }
    }
    // set up cplex solver
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.solve();
    // obtain optimial solution of SAA problem
    std::vector<double> x_est(x_size,0.0);
    for (int x_index = 0; x_index < x_size; ++x_index) {
        x_est[x_index] = cplex.getValue(x[x_index]);
    }
    return x_est;
}

validationResult oneStagePiecewiseLP_SAA_validation(const std::string& folder_path, const std::vector<double>& x_candidate) {
    // STEP 1: INITIALIZATION
    // create directory paths for database and model
    std::string DB_path = folder_path + "/DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(SAA).txt";
    // create database
    std::vector<std::vector<dataPoint>> rawObjectiveDB = readNonparametricDB(DB_path);
    // create model structure
    oneStageParameters parameters = readOneStageParameters(model_path);
    // number of scenarios
    long num_scenarios = rawObjectiveDB[0].size();
    // size of x
    long x_size = parameters.x.size();
    // create sto object
    oneStagePiecewiseLinearObjective parameters_objective = readStochasticOneStagePiecewiseLP(sto_path);
    // merge database and parameters of objeective
    oneStagePiecewiseLinearObjectiveDB objectiveDB = mergeDB_randomObjective(rawObjectiveDB, parameters_objective, x_size);
    // calculate cost in each scenario
    double mean_cost = 0;
    double meanOfSquare_cost = 0;
    for (int scenario_index = 0; scenario_index < num_scenarios; ++scenario_index) {
        double tempCost = 0; // cost for this scenario
        // c1^T X
        for (auto x_iterator : objectiveDB.c1[0][scenario_index]) {
            int x_index = x_iterator.first; // index of x
            double coefficient = x_iterator.second; // coeffcient for x[x_index]
            tempCost += coefficient * x_candidate[x_index];
        }
        // e1
        tempCost += objectiveDB.e1[0][scenario_index];
        // PL
        double PLCost1 = 0;
        // c2^T x
        for (auto x_iterator : objectiveDB.c2[0][scenario_index]) {
            int x_index = x_iterator.first; // index of x
            double coefficient = x_iterator.second; // coeffcient for x[x_index]
            PLCost1 += coefficient * x_candidate[x_index];
        }
        PLCost1 += objectiveDB.e2[0][scenario_index];
        double PLCost2 = 0;
        // c2^T x
        for (auto x_iterator : objectiveDB.c3[0][scenario_index]) {
            int x_index = x_iterator.first; // index of x
            double coefficient = x_iterator.second; // coeffcient for x[x_index]
            PLCost2 += coefficient * x_candidate[x_index];
        }
        PLCost2 += objectiveDB.e3[0][scenario_index];
        if (PLCost1 > PLCost2) {
            tempCost += objectiveDB.cPiecewise * PLCost1;
        }
        else {
            tempCost += objectiveDB.cPiecewise * PLCost2;
        }
        // finish calculating tempCost
        mean_cost += tempCost / num_scenarios;
        meanOfSquare_cost = (meanOfSquare_cost * scenario_index + tempCost * tempCost) / ((double) scenario_index + 1);
    }
    double variance_cost = meanOfSquare_cost * (num_scenarios / ((double)num_scenarios - 1)) - mean_cost * (num_scenarios / ((double) num_scenarios - 1));
    validationResult result;
    result.mean = mean_cost;
    result.variance = variance_cost;
    double halfMargin = result.Zalpha * sqrt(result.variance / num_scenarios);
    result.CI_lower = mean_cost - halfMargin;
    result.CI_upper = mean_cost + halfMargin;
    std::cout << "Mean of cost: " << mean_cost << std::endl;
    std::cout << "Variance of cost: " << variance_cost << std::endl;
    std::cout << result.alpha << "% CI Lower Bound: " << result.CI_lower << std::endl;
    std::cout << result.alpha << "% CI Upper Bound: " << result.CI_upper << std::endl;
    return result;
}
