//
//  twoStageLP_NSQG.cpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/23/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#include "twoStageLP_NSQG.hpp"

// random quantities only appear on the right hand side of constraints in the second stage problem
std::vector<double> twoStageLP_random_b(const std::string& folder_path, int max_iterations, double initial_stepsize, std::vector<double> x_init) {
    // STEP 1: INITIALIZATION
    // create directory paths for database and model
    std::string be_DB_path = folder_path + "/be_DB.txt";
    std::string bi_DB_path = folder_path + "/bi_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    // convert all the paths into constant chars
    const char* be_DB_path_const = be_DB_path.c_str();
    const char* bi_DB_path_const = bi_DB_path.c_str();
    const char* model_path_const = model_path.c_str();
    const char* sto_path_const = sto_path.c_str();
    // create stream object
    std::ifstream readFile_be(be_DB_path_const);
    std::ifstream readFile_bi(bi_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> be_DB;
    std::vector<std::vector<dataPoint>> bi_DB;
    // create model structure
    twoStageParameters model_parameters;
    // create sto object
    secondStageRHS bRHS;
    // read  be
    if (readFile_be.is_open()) {
        std::cout << "be_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read be database
        be_DB = readNonparametricDB(be_DB_path);
    }
    else {
        readFile_be.close(); // close the file
        std::cout << "be_DB is not found!" << std::endl;
    }
    // read bi
    if (readFile_bi.is_open()) {
        std::cout << "bi_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read bi database
        bi_DB = readNonparametricDB(bi_DB_path);
    }
    else {
        readFile_bi.close(); // close the file
        std::cout << "bi_DB is not found!" << std::endl;
    }
    // read model
    model_parameters = readTwoStageParameters(model_path);
    // read sto object
    bRHS = readStochastic(sto_path);
    // create be and bi in std::vector<double> form
    secondStageRHSDB bRHSDB = mergeDB_randomVector(be_DB, bi_DB, bRHS); //
    /*
    // for debugging
    std::cout << ">>>>>bi_DB<<<<<" << std::endl;
    printNonparametricDB(bi_DB);
    std::cout << ">>>>>bi_DB<<<<<" << std::endl;
    std::cout << ">>>>>Model<<<<<" << std::endl;
    printTwoStageParameters(model_parameters);
    std::cout << ">>>>>Model<<<<<" << std::endl;
    std::cout << ">>>>>Sto<<<<<" << std::endl;
    printStochastic(bRHS);
    std::cout << ">>>>>Sto<<<<<" << std::endl;
     */
    std::cout << "STEP 2: SOLVE TWO STAGE LINEAR PROGRAMMING" << std::endl;
    // STEP 2: SOLVE TWO STAGE LINEAR PROGRAMMING
    std::vector<double> x_est;
    x_est = twoStageLP_solver(model_parameters, bRHSDB.be_database, bRHSDB.bi_database, bRHSDB.weight_database, initial_stepsize, max_iterations, x_init);
    // STEP 3: OUTPUT RESULT
    return x_est;
}

// computational results will be stored in text file
std::vector<double> twoStageLP_random_b_outputResults(const std::string& folder_path, int max_iterations, double initial_stepsize, std::vector<double> x_init) {
    // STEP 1: INITIALIZATION
    // create directory paths for database and model
    std::string be_DB_path = folder_path + "/be_DB.txt";
    std::string bi_DB_path = folder_path + "/bi_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults.txt";
    // convert all the paths into constant chars
    const char* be_DB_path_const = be_DB_path.c_str();
    const char* bi_DB_path_const = bi_DB_path.c_str();
    const char* model_path_const = model_path.c_str();
    const char* sto_path_const = sto_path.c_str();
    const char* resultsOutput_path_const = resultsOutput_path.c_str();
    // create stream object
    std::ifstream readFile_be(be_DB_path_const);
    std::ifstream readFile_bi(bi_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> be_DB;
    std::vector<std::vector<dataPoint>> bi_DB;
    // create model structure
    twoStageParameters model_parameters;
    // create sto object
    secondStageRHS bRHS;
    // read  be
    if (readFile_be.is_open()) {
        std::cout << "be_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read be database
        be_DB = readNonparametricDB(be_DB_path);
    }
    else {
        readFile_be.close(); // close the file
        std::cout << "be_DB is not found!" << std::endl;
    }
    // read bi
    if (readFile_bi.is_open()) {
        std::cout << "bi_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read bi database
        bi_DB = readNonparametricDB(bi_DB_path);
    }
    else {
        readFile_bi.close(); // close the file
        std::cout << "bi_DB is not found!" << std::endl;
    }
    // read model
    model_parameters = readTwoStageParameters(model_path);
    // read sto object
    bRHS = readStochastic(sto_path);
    // create be and bi in std::vector<double> form
    secondStageRHSDB bRHSDB = mergeDB_randomVector(be_DB, bi_DB, bRHS); //
    // for debugging
    /*
    std::cout << ">>>>>bi_DB<<<<<" << std::endl;
    printNonparametricDB(bi_DB);
    std::cout << ">>>>>bi_DB<<<<<" << std::endl;
    std::cout << ">>>>>Model<<<<<" << std::endl;
    printTwoStageParameters(model_parameters);
    std::cout << ">>>>>Model<<<<<" << std::endl;
    std::cout << ">>>>>Sto<<<<<" << std::endl;
    printStochastic(bRHS);
    std::cout << ">>>>>Sto<<<<<" << std::endl;
     */
    std::cout << "STEP 2: SOLVE TWO STAGE LINEAR PROGRAMMING" << std::endl;
    // STEP 2: SOLVE TWO STAGE LINEAR PROGRAMMING
    std::vector<double> x_est;
    x_est = twoStageLP_solver(model_parameters, bRHSDB.be_database, bRHSDB.bi_database, bRHSDB.weight_database, initial_stepsize, max_iterations, x_init, resultsOutput_path_const);
    // STEP 3: OUTPUT RESULT
    return x_est;
}

// Robust Nonparametric SQG
std::vector<double> robustTwoStageLP_random_b_outputResults(const std::string& folder_path, int max_iterations, std::vector<double> x_init, double Dx, int m, double M) {
    // m multiplier for iterates
    // M upper bound of the L2 norm of quasi-gradient
    // STEP 1: INITIALIZATION
    // create directory paths for database and model
    std::string be_DB_path = folder_path + "/be_DB.txt";
    std::string bi_DB_path = folder_path + "/bi_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    std::string resultsOutput_path = folder_path + "/computationalResults(Robust_NSQG).txt";
    // convert all the paths into constant chars
    const char* be_DB_path_const = be_DB_path.c_str();
    const char* bi_DB_path_const = bi_DB_path.c_str();
    const char* model_path_const = model_path.c_str();
    const char* sto_path_const = sto_path.c_str();
    const char* resultsOutput_path_const = resultsOutput_path.c_str();
    // create stream object
    std::ifstream readFile_be(be_DB_path_const);
    std::ifstream readFile_bi(bi_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> be_DB;
    std::vector<std::vector<dataPoint>> bi_DB;
    // create model structure
    twoStageParameters model_parameters;
    // create sto object
    secondStageRHS bRHS;
    // read  be
    if (readFile_be.is_open()) {
        std::cout << "be_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read be database
        be_DB = readNonparametricDB(be_DB_path);
    }
    else {
        readFile_be.close(); // close the file
        std::cout << "be_DB is not found!" << std::endl;
    }
    // read bi
    if (readFile_bi.is_open()) {
        std::cout << "bi_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read bi database
        bi_DB = readNonparametricDB(bi_DB_path);
    }
    else {
        readFile_bi.close(); // close the file
        std::cout << "bi_DB is not found!" << std::endl;
    }
    // read model
    model_parameters = readTwoStageParameters(model_path);
    // read sto object
    bRHS = readStochastic(sto_path);
    // create be and bi in std::vector<double> form
    secondStageRHSDB bRHSDB = mergeDB_randomVector(be_DB, bi_DB, bRHS); //
    // for debugging
    /*
    std::cout << ">>>>>bi_DB<<<<<" << std::endl;
    printNonparametricDB(bi_DB);
    std::cout << ">>>>>bi_DB<<<<<" << std::endl;
    std::cout << ">>>>>Model<<<<<" << std::endl;
    printTwoStageParameters(model_parameters);
    std::cout << ">>>>>Model<<<<<" << std::endl;
    std::cout << ">>>>>Sto<<<<<" << std::endl;
    printStochastic(bRHS);
    std::cout << ">>>>>Sto<<<<<" << std::endl;
     */
    std::cout << "STEP 2: SOLVE TWO STAGE LINEAR PROGRAMMING" << std::endl;
    // STEP 2: SOLVE TWO STAGE LINEAR PROGRAMMING
    std::vector<double> x_est;
    x_est = robustTwoStageLP_solver(model_parameters, bRHSDB.be_database, bRHSDB.bi_database, bRHSDB.weight_database, max_iterations, x_init, Dx, m, M, resultsOutput_path_const);
    // STEP 3: OUTPUT RESULT
    return x_est;
}

// estimate cost via validation set
double twoStageLP_validation(const std::string& folder_path, const std::vector<double>& x_candidate) {
    // STEP 1: INITIALIZATION
    // create directory paths for database and model
    std::string be_DB_path = folder_path + "/be_DB.txt";
    std::string bi_DB_path = folder_path + "/bi_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    // convert all the paths into constant chars
    const char* be_DB_path_const = be_DB_path.c_str();
    const char* bi_DB_path_const = bi_DB_path.c_str();
    const char* model_path_const = model_path.c_str();
    const char* sto_path_const = sto_path.c_str();
    // create stream object
    std::ifstream readFile_be(be_DB_path_const);
    std::ifstream readFile_bi(bi_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> be_DB;
    std::vector<std::vector<dataPoint>> bi_DB;
    // create model structure
    twoStageParameters model_parameters;
    // create sto object
    secondStageRHS bRHS;
    // read  be
    if (readFile_be.is_open()) {
        std::cout << "be_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read be database
        be_DB = readNonparametricDB(be_DB_path);
    }
    else {
        readFile_be.close(); // close the file
        std::cout << "be_DB is not found!" << std::endl;
    }
    // read bi
    if (readFile_bi.is_open()) {
        std::cout << "bi_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read bi database
        bi_DB = readNonparametricDB(bi_DB_path);
    }
    else {
        readFile_bi.close(); // close the file
        std::cout << "bi_DB is not found!" << std::endl;
    }
    // read model
    model_parameters = readTwoStageParameters(model_path);
    // read sto object
    bRHS = readStochastic(sto_path);
    // create be and bi in std::vector<double> form
    secondStageRHSDB bRHSDB = mergeDB_randomVector(be_DB, bi_DB, bRHS); //
    // for debugging
    /*
     std::cout << ">>>>>bi_DB<<<<<" << std::endl;
     printNonparametricDB(bi_DB);
     std::cout << ">>>>>bi_DB<<<<<" << std::endl;
     std::cout << ">>>>>Model<<<<<" << std::endl;
     printTwoStageParameters(model_parameters);
     std::cout << ">>>>>Model<<<<<" << std::endl;
     std::cout << ">>>>>Sto<<<<<" << std::endl;
     printStochastic(bRHS);
     std::cout << ">>>>>Sto<<<<<" << std::endl;
     */
    // size of the validation set
    // size of be
    long be_database_size = bRHSDB.be_database.size();
    long be_dataset_size = 0;
    if (be_database_size > 0) {
        be_dataset_size = bRHSDB.be_database[0].size();
    }
    // size of bi
    long bi_database_size = bRHSDB.bi_database.size();
    long bi_dataset_size = 0;
    if (bi_database_size > 0) {
        bi_dataset_size = bRHSDB.bi_database[0].size();
    }
    // determine appropriate model
    double secondStageTotalCost = 0;
    if (be_dataset_size > 0 && bi_dataset_size > 0) { // second stage subproblem has both equality constraints and inequality constraints
        for (int scenario_index = 0; scenario_index < bi_dataset_size; ++scenario_index) {
            double tempCost = twoStageLP_secondStageCost(x_candidate, model_parameters.d, model_parameters.De, model_parameters.Ce,bRHSDB.be_database[0][scenario_index], model_parameters.Di, model_parameters.Ci, bRHSDB.bi_database[0][scenario_index]);
            secondStageTotalCost += tempCost;
        }
        secondStageTotalCost = secondStageTotalCost / bi_dataset_size;
    }
    else if (be_dataset_size > 0 && bi_dataset_size < 1) { // second stage subproblem only has equality constraints
        for (int scenario_index = 0; scenario_index < be_dataset_size; ++scenario_index) {
            double tempCost = twoStageLP_secondStageCost_equality(x_candidate, model_parameters.d, model_parameters.De, model_parameters.Ce,bRHSDB.be_database[0][scenario_index]);
            secondStageTotalCost += tempCost;
        }
        secondStageTotalCost = secondStageTotalCost / be_dataset_size;
    }
    else if (be_dataset_size < 1 && bi_dataset_size > 0) { // second stage subproblem only has inequality constraints
        for (int scenario_index = 0; scenario_index < bi_dataset_size; ++scenario_index) {
            double tempCost = twoStageLP_secondStageCost_inequality(x_candidate, model_parameters.d, model_parameters.Di, model_parameters.Ci,bRHSDB.bi_database[0][scenario_index]);
            secondStageTotalCost += tempCost;
        }
        secondStageTotalCost = secondStageTotalCost / bi_dataset_size;
    }
    else { // second stage subproblem is unconstrained
        std::cout << "Warning: Second stage subproblem is unconstrained! Return 0." << std::endl;
    }
    double firstStageCost = model_parameters.c * x_candidate;
    return firstStageCost + secondStageTotalCost;
}

double twoStageLP_validation_outputResults(const std::string& folder_path, const std::vector<double>& x_candidate) {
    // STEP 1: INITIALIZATION
    // create directory paths for database and model
    std::string be_DB_path = folder_path + "/be_DB.txt";
    std::string bi_DB_path = folder_path + "/bi_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    // convert all the paths into constant chars
    const char* be_DB_path_const = be_DB_path.c_str();
    const char* bi_DB_path_const = bi_DB_path.c_str();
    const char* model_path_const = model_path.c_str();
    const char* sto_path_const = sto_path.c_str();
    // create stream object
    std::ifstream readFile_be(be_DB_path_const);
    std::ifstream readFile_bi(bi_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> be_DB;
    std::vector<std::vector<dataPoint>> bi_DB;
    // create model structure
    twoStageParameters model_parameters;
    // create sto object
    secondStageRHS bRHS;
    // read  be
    if (readFile_be.is_open()) {
        std::cout << "be_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read be database
        be_DB = readNonparametricDB(be_DB_path);
    }
    else {
        readFile_be.close(); // close the file
        std::cout << "be_DB is not found!" << std::endl;
    }
    // read bi
    if (readFile_bi.is_open()) {
        std::cout << "bi_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read bi database
        bi_DB = readNonparametricDB(bi_DB_path);
    }
    else {
        readFile_bi.close(); // close the file
        std::cout << "bi_DB is not found!" << std::endl;
    }
    // read model
    model_parameters = readTwoStageParameters(model_path);
    // read sto object
    bRHS = readStochastic(sto_path);
    // create be and bi in std::vector<double> form
    secondStageRHSDB bRHSDB = mergeDB_randomVector(be_DB, bi_DB, bRHS); //
    // for debugging
    /*
     std::cout << ">>>>>bi_DB<<<<<" << std::endl;
     printNonparametricDB(bi_DB);
     std::cout << ">>>>>bi_DB<<<<<" << std::endl;
     std::cout << ">>>>>Model<<<<<" << std::endl;
     printTwoStageParameters(model_parameters);
     std::cout << ">>>>>Model<<<<<" << std::endl;
     std::cout << ">>>>>Sto<<<<<" << std::endl;
     printStochastic(bRHS);
     std::cout << ">>>>>Sto<<<<<" << std::endl;
     */
    // size of the validation set
    // size of be
    long be_database_size = bRHSDB.be_database.size();
    long be_dataset_size = 0;
    if (be_database_size > 0) {
        be_dataset_size = bRHSDB.be_database[0].size();
    }
    // size of bi
    long bi_database_size = bRHSDB.bi_database.size();
    long bi_dataset_size = 0;
    if (bi_database_size > 0) {
        bi_dataset_size = bRHSDB.bi_database[0].size();
    }
    // determine appropriate model
    double secondStageTotalCost = 0;
    if (be_dataset_size > 0 && bi_dataset_size > 0) { // second stage subproblem has both equality constraints and inequality constraints
        for (int scenario_index = 0; scenario_index < bi_dataset_size; ++scenario_index) {
            double tempCost = twoStageLP_secondStageCost(x_candidate, model_parameters.d, model_parameters.De, model_parameters.Ce,bRHSDB.be_database[0][scenario_index], model_parameters.Di, model_parameters.Ci, bRHSDB.bi_database[0][scenario_index]);
            secondStageTotalCost += tempCost / bi_dataset_size;
        }
    }
    else if (be_dataset_size > 0 && bi_dataset_size < 1) { // second stage subproblem only has equality constraints
        for (int scenario_index = 0; scenario_index < be_dataset_size; ++scenario_index) {
            double tempCost = twoStageLP_secondStageCost_equality(x_candidate, model_parameters.d, model_parameters.De, model_parameters.Ce,bRHSDB.be_database[0][scenario_index]);
            secondStageTotalCost += tempCost / be_dataset_size;
        }
    }
    else if (be_dataset_size < 1 && bi_dataset_size > 0) { // second stage subproblem only has inequality constraints
        for (int scenario_index = 0; scenario_index < bi_dataset_size; ++scenario_index) {
            double tempCost = twoStageLP_secondStageCost_inequality(x_candidate, model_parameters.d, model_parameters.Di, model_parameters.Ci,bRHSDB.bi_database[0][scenario_index]);
            secondStageTotalCost += tempCost / bi_dataset_size;
        }
    }
    else { // second stage subproblem is unconstrained
        std::cout << "Warning: Second stage subproblem is unconstrained! Return 0." << std::endl;
    }
    double firstStageCost = model_parameters.c * x_candidate;
    // write computational results
    std::string outputResults_path = folder_path + "/validationResults.txt";
    const char* outputResults_path_const = outputResults_path.c_str();
    std::fstream writeFile;
    writeFile.open(outputResults_path_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "***************************************************\n";
    writeFile << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution: ";
    long x_size = x_candidate.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_candidate[index] << ", ";
    }
    writeFile << x_candidate[x_size - 1] << "\n";
    writeFile << "Number of data points: " << max(bi_dataset_size,be_dataset_size) << "\n";
    writeFile << "Validation cost: " << firstStageCost + secondStageTotalCost << "\n";
    writeFile << "***************************************************\n";
    return firstStageCost + secondStageTotalCost;
}

// get the cost of second stage subproblem
double twoStageLP_secondStageCost(const std::vector<double>& x, const std::vector<double>& d, const std::vector<std::vector<double>>& De, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be, const std::vector<std::vector<double>>& Di, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi) {
    // derive the sizes of input parameters
    long d_size = d.size();
    long De_rowsize = De.size();
    long De_colsize = 0;
    if (De_rowsize > 0) {
        De_colsize = De[0].size();
    }
    long Ce_rowsize = Ce.size();
    long Ce_colsize = 0;
    if (Ce_rowsize > 0) {
        Ce_colsize = Ce[0].size();
    }
    long be_rowsize = be.size();
    long Di_rowsize = Di.size();
    long Di_colsize = 0;
    if (Di_rowsize > 0) {
        Di_colsize = Di[0].size();
    }
    long Ci_rowsize = Ci.size();
    long Ci_colsize = 0;
    if (Ci_rowsize > 0) {
        Ci_colsize = Ci[0].size();
    }
    long bi_size = bi.size();
    long equality_size = De_rowsize; // number of equality constraints
    long inequality_size = Di_rowsize; // number of inequality constraints
    // set up the model
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray y(env,d_size,-IloInfinity,IloInfinity,ILOFLOAT); // second stage decision variables
    mod.add(y);
    // objective function
    IloExpr expr_obj(env);
    for (int d_index = 0; d_index < d_size; ++d_index) {
        expr_obj += d[d_index] * y[d_index];
    }
    IloObjective obj = IloMinimize(env,expr_obj);
    mod.add(obj);
    // constraints
    IloRangeArray constraintsEquality(env);
    // equality constraints
    for (int equality_index = 0; equality_index < equality_size; ++equality_index) {
        IloExpr expr_equality(env);
        for (int De_index = 0; De_index < De_colsize; ++De_index) {
            expr_equality += De[equality_index][De_index] * y[De_index];
        }
        double Ce_times_x = Ce[equality_index] * x;
        expr_equality += Ce_times_x - be[equality_index];
        constraintsEquality.add(expr_equality == 0);
    }
    if (equality_size > 0) { // if there is at least one equality constraint
        mod.add(constraintsEquality);
    }
    // inequality constraints
    IloRangeArray constraintsInequality(env);
    for (int inequality_index = 0; inequality_index < inequality_size; ++inequality_index) {
        IloExpr expr_inequality(env);
        for (int Di_index = 0; Di_index < Di_colsize; ++Di_index) {
            expr_inequality += Di[inequality_index][Di_index] * y[Di_index];
        }
        double Ci_times_x = Ci[inequality_index] * x;
        expr_inequality += Ci_times_x - bi[inequality_index];
        constraintsInequality.add(expr_inequality <= 0);
    }
    if (inequality_size > 0) { // if there is at least one inequality constraint
        mod.add(constraintsInequality);
    }
    // set up cplex solver
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setOut(env.getNullStream());
    IloBool solvable_flag = cplex.solve();
    IloNumArray dual_equality(env);
    IloNumArray dual_inequality(env);
    double secondStageCost = cplex.getObjValue();
    env.end(); // end the environment
    return secondStageCost;
}


double twoStageLP_secondStageCost_inequality(const std::vector<double>& x, const std::vector<double>& d, const std::vector<std::vector<double>>& Di, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi) {
    // derive the sizes of input parameters
    long d_size = d.size();
    long Di_rowsize = Di.size();
    long Di_colsize = 0;
    if (Di_rowsize > 0) {
        Di_colsize = Di[0].size();
    }
    long Ci_rowsize = Ci.size();
    long Ci_colsize = 0;
    if (Ci_rowsize > 0) {
        Ci_colsize = Ci[0].size();
    }
    long bi_size = bi.size();
    long inequality_size = Di_rowsize; // number of inequality constraints
    // set up the model
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray y(env,d_size,-IloInfinity,IloInfinity,ILOFLOAT); // second stage decision variables
    mod.add(y);
    // objective function
    IloExpr expr_obj(env);
    for (int d_index = 0; d_index < d_size; ++d_index) {
        expr_obj += d[d_index] * y[d_index];
    }
    IloObjective obj = IloMinimize(env,expr_obj);
    mod.add(obj);
    // constraints
    // inequality constraints
    IloRangeArray constraintsInequality(env);
    for (int inequality_index = 0; inequality_index < inequality_size; ++inequality_index) {
        IloExpr expr_inequality(env);
        for (int Di_index = 0; Di_index < Di_colsize; ++Di_index) {
            expr_inequality += Di[inequality_index][Di_index] * y[Di_index];
        }
        double Ci_times_x = Ci[inequality_index] * x;
        expr_inequality += Ci_times_x - bi[inequality_index];
        constraintsInequality.add(expr_inequality <= 0);
    }
    if (inequality_size > 0) { // if there is at least one inequality constraint
        mod.add(constraintsInequality);
    }
    // set up cplex solver
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setOut(env.getNullStream());
    IloBool solvable_flag = cplex.solve();
    IloNumArray dual_equality(env);
    IloNumArray dual_inequality(env);
    double secondStageCost = cplex.getObjValue();
    env.end(); // end the environment
    return secondStageCost;
}


double twoStageLP_secondStageCost_equality(const std::vector<double>& x, const std::vector<double>& d, const std::vector<std::vector<double>>& De, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be) {
    // derive the sizes of input parameters
    long d_size = d.size();
    long De_rowsize = De.size();
    long De_colsize = 0;
    if (De_rowsize > 0) {
        De_colsize = De[0].size();
    }
    long Ce_rowsize = Ce.size();
    long Ce_colsize = 0;
    if (Ce_rowsize > 0) {
        Ce_colsize = Ce[0].size();
    }
    long be_rowsize = be.size();
    long equality_size = De_rowsize; // number of equality constraints
    // set up the model
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray y(env,d_size,-IloInfinity,IloInfinity,ILOFLOAT); // second stage decision variables
    mod.add(y);
    // objective function
    IloExpr expr_obj(env);
    for (int d_index = 0; d_index < d_size; ++d_index) {
        expr_obj += d[d_index] * y[d_index];
    }
    IloObjective obj = IloMinimize(env,expr_obj);
    mod.add(obj);
    // constraints
    IloRangeArray constraintsEquality(env);
    // equality constraints
    for (int equality_index = 0; equality_index < equality_size; ++equality_index) {
        IloExpr expr_equality(env);
        for (int De_index = 0; De_index < De_colsize; ++De_index) {
            expr_equality += De[equality_index][De_index] * y[De_index];
        }
        double Ce_times_x = Ce[equality_index] * x;
        expr_equality += Ce_times_x - be[equality_index];
        constraintsEquality.add(expr_equality == 0);
    }
    if (equality_size > 0) { // if there is at least one equality constraint
        mod.add(constraintsEquality);
    }
    // set up cplex solver
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setOut(env.getNullStream());
    IloBool solvable_flag = cplex.solve();
    IloNumArray dual_equality(env);
    IloNumArray dual_inequality(env);
    double secondStageCost = cplex.getObjValue();
    env.end(); // end the environment
    return secondStageCost;
}

// estimate the value of nonparametric estimator
double twoStageLP_estimator(const std::string& folder_path, const std::vector<double>& x_candidate, int dataset_index) {
    // STEP 1: INITIALIZATION
    // create directory paths for database and model
    std::string be_DB_path = folder_path + "/be_DB.txt";
    std::string bi_DB_path = folder_path + "/bi_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    // convert all the paths into constant chars
    const char* be_DB_path_const = be_DB_path.c_str();
    const char* bi_DB_path_const = bi_DB_path.c_str();
    const char* model_path_const = model_path.c_str();
    const char* sto_path_const = sto_path.c_str();
    // create stream object
    std::ifstream readFile_be(be_DB_path_const);
    std::ifstream readFile_bi(bi_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> be_DB;
    std::vector<std::vector<dataPoint>> bi_DB;
    // create model structure
    twoStageParameters model_parameters;
    // create sto object
    secondStageRHS bRHS;
    // read  be
    std::cout << be_DB_path_const << std::endl;
    std::cout << bi_DB_path_const << std::endl;
    if (readFile_be.is_open()) {
        std::cout << "be_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read be database
        be_DB = readNonparametricDB(be_DB_path);
    }
    else {
        readFile_be.close(); // close the file
        std::cout << "be_DB is not found!" << std::endl;
    }
    // read bi
    if (readFile_bi.is_open()) {
        std::cout << "bi_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read bi database
        bi_DB = readNonparametricDB(bi_DB_path);
    }
    else {
        readFile_bi.close(); // close the file
        std::cout << "bi_DB is not found!" << std::endl;
    }
    // read model
    model_parameters = readTwoStageParameters(model_path);
    // read sto object
    bRHS = readStochastic(sto_path);
    // create be and bi in std::vector<double> form
    secondStageRHSDB bRHSDB = mergeDB_randomVector(be_DB, bi_DB, bRHS); //
    // for debugging
    /*
     std::cout << ">>>>>bi_DB<<<<<" << std::endl;
     printNonparametricDB(bi_DB);
     std::cout << ">>>>>bi_DB<<<<<" << std::endl;
     std::cout << ">>>>>Model<<<<<" << std::endl;
     printTwoStageParameters(model_parameters);
     std::cout << ">>>>>Model<<<<<" << std::endl;
     std::cout << ">>>>>Sto<<<<<" << std::endl;
     printStochastic(bRHS);
     std::cout << ">>>>>Sto<<<<<" << std::endl;
     */
    // size of the validation set
    // size of be
    long be_database_size = bRHSDB.be_database.size();
    long be_dataset_size = 0;
    if (be_database_size > 0) {
        be_dataset_size = bRHSDB.be_database[dataset_index].size();
    }
    // size of bi
    long bi_database_size = bRHSDB.bi_database.size();
    long bi_dataset_size = 0;
    if (bi_database_size > 0) {
        bi_dataset_size = bRHSDB.bi_database[dataset_index].size();
    }
    // determine appropriate model
    double secondStageTotalCost = 0;
    double firstStageCost = model_parameters.c * x_candidate;
    double mean = 0;
    int num_dataPoints = 0;
    if (be_dataset_size > 0 && bi_dataset_size > 0) { // second stage subproblem has both equality constraints and inequality constraints
        num_dataPoints = be_dataset_size;
        for (int scenario_index = 0; scenario_index < bi_dataset_size; ++scenario_index) {
            double tempCost = twoStageLP_secondStageCost(x_candidate, model_parameters.d, model_parameters.De, model_parameters.Ce,bRHSDB.be_database[dataset_index][scenario_index], model_parameters.Di, model_parameters.Ci, bRHSDB.bi_database[dataset_index][scenario_index]);
            secondStageTotalCost += tempCost * bRHSDB.weight_database[dataset_index][scenario_index];
        }
    }
    else if (be_dataset_size > 0 && bi_dataset_size < 1) { // second stage subproblem only has equality constraints
        num_dataPoints = be_dataset_size;
        for (int scenario_index = 0; scenario_index < be_dataset_size; ++scenario_index) {
            double tempCost = twoStageLP_secondStageCost_equality(x_candidate, model_parameters.d, model_parameters.De, model_parameters.Ce,bRHSDB.be_database[dataset_index][scenario_index]);
            secondStageTotalCost += tempCost * bRHSDB.weight_database[dataset_index][scenario_index];
        }
    }
    else if (be_dataset_size < 1 && bi_dataset_size > 0) { // second stage subproblem only has inequality constraints
        num_dataPoints = bi_dataset_size;
        for (int scenario_index = 0; scenario_index < bi_dataset_size; ++scenario_index) {
            double tempCost = twoStageLP_secondStageCost_inequality(x_candidate, model_parameters.d, model_parameters.Di, model_parameters.Ci,bRHSDB.bi_database[dataset_index][scenario_index]);
            secondStageTotalCost += tempCost * bRHSDB.weight_database[dataset_index][scenario_index];
        }
    }
    else { // second stage subproblem is unconstrained
        std::cout << "Warning: Second stage subproblem is unconstrained! Return 0." << std::endl;
        return 0;
    }
    // write computational results
    std::string outputResults_path = folder_path + "/estimatorResults.txt";
    const char* outputResults_path_const = outputResults_path.c_str();
    std::fstream writeFile;
    writeFile.open(outputResults_path_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "***************************************************\n";
    writeFile << "Folder Path                    : " << folder_path << std::endl;
    writeFile << "Estimating the value of estimator at the the point below\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution             : ";
    long x_size = x_candidate.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_candidate[index] << ", ";
    }
    writeFile << x_candidate[x_size - 1] << "\n";
    writeFile << "Number of Scenarios            : " << num_dataPoints << std::endl;
    writeFile << "Dataset Index                  : " << dataset_index << std::endl;
    writeFile << "================Value of Estimated Function==========\n";
    writeFile << firstStageCost + secondStageTotalCost << std::endl;
    writeFile << "***************************************************\n";
    double result = firstStageCost + secondStageTotalCost;
    return result;
}

// estimate the quality of candiate solution (v2)
validationResult twoStageLP_validation_outputResultsV2(const std::string& folder_path, const std::vector<double>& x_candidate) {
    // STEP 1: INITIALIZATION
    // create directory paths for database and model
    std::string be_DB_path = folder_path + "/be_DB.txt";
    std::string bi_DB_path = folder_path + "/bi_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    // convert all the paths into constant chars
    const char* be_DB_path_const = be_DB_path.c_str();
    const char* bi_DB_path_const = bi_DB_path.c_str();
    const char* model_path_const = model_path.c_str();
    const char* sto_path_const = sto_path.c_str();
    // create stream object
    std::ifstream readFile_be(be_DB_path_const);
    std::ifstream readFile_bi(bi_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> be_DB;
    std::vector<std::vector<dataPoint>> bi_DB;
    // create model structure
    twoStageParameters model_parameters;
    // create sto object
    secondStageRHS bRHS;
    // read  be
    std::cout << be_DB_path_const << std::endl;
    std::cout << bi_DB_path_const << std::endl;
    if (readFile_be.is_open()) {
        std::cout << "be_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read be database
        be_DB = readNonparametricDB(be_DB_path);
    }
    else {
        readFile_be.close(); // close the file
        std::cout << "be_DB is not found!" << std::endl;
    }
    // read bi
    if (readFile_bi.is_open()) {
        std::cout << "bi_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read bi database
        bi_DB = readNonparametricDB(bi_DB_path);
    }
    else {
        readFile_bi.close(); // close the file
        std::cout << "bi_DB is not found!" << std::endl;
    }
    // read model
    model_parameters = readTwoStageParameters(model_path);
    // read sto object
    bRHS = readStochastic(sto_path);
    // create be and bi in std::vector<double> form
    secondStageRHSDB bRHSDB = mergeDB_randomVector(be_DB, bi_DB, bRHS); //
    // for debugging
    /*
     std::cout << ">>>>>bi_DB<<<<<" << std::endl;
     printNonparametricDB(bi_DB);
     std::cout << ">>>>>bi_DB<<<<<" << std::endl;
     std::cout << ">>>>>Model<<<<<" << std::endl;
     printTwoStageParameters(model_parameters);
     std::cout << ">>>>>Model<<<<<" << std::endl;
     std::cout << ">>>>>Sto<<<<<" << std::endl;
     printStochastic(bRHS);
     std::cout << ">>>>>Sto<<<<<" << std::endl;
     */
    // size of the validation set
    // size of be
    long be_database_size = bRHSDB.be_database.size();
    long be_dataset_size = 0;
    if (be_database_size > 0) {
        be_dataset_size = bRHSDB.be_database[0].size();
    }
    // size of bi
    long bi_database_size = bRHSDB.bi_database.size();
    long bi_dataset_size = 0;
    if (bi_database_size > 0) {
        bi_dataset_size = bRHSDB.bi_database[0].size();
    }
    // determine appropriate model
    double secondStageTotalCost = 0;
    double firstStageCost = model_parameters.c * x_candidate;
    double variance_P = 0; // intermediate component for calculating variance
    if (be_dataset_size > 0 && bi_dataset_size > 0) { // second stage subproblem has both equality constraints and inequality constraints
        for (int scenario_index = 0; scenario_index < bi_dataset_size; ++scenario_index) {
            double tempCost = twoStageLP_secondStageCost(x_candidate, model_parameters.d, model_parameters.De, model_parameters.Ce,bRHSDB.be_database[0][scenario_index], model_parameters.Di, model_parameters.Ci, bRHSDB.bi_database[0][scenario_index]);
            secondStageTotalCost += tempCost / bi_dataset_size;
            double tempTotalCost = tempCost + firstStageCost;
            double _n_ = scenario_index + 1;
            variance_P = variance_P * (_n_ - 1) / _n_ + tempTotalCost * tempTotalCost / _n_;
        }
    }
    else if (be_dataset_size > 0 && bi_dataset_size < 1) { // second stage subproblem only has equality constraints
        for (int scenario_index = 0; scenario_index < be_dataset_size; ++scenario_index) {
            double tempCost = twoStageLP_secondStageCost_equality(x_candidate, model_parameters.d, model_parameters.De, model_parameters.Ce,bRHSDB.be_database[0][scenario_index]);
            secondStageTotalCost += tempCost / be_dataset_size;
            double tempTotalCost = tempCost + firstStageCost;
            double _n_ = scenario_index + 1;
            variance_P = variance_P * (_n_ - 1) / _n_ + tempTotalCost * tempTotalCost / _n_;
        }
    }
    else if (be_dataset_size < 1 && bi_dataset_size > 0) { // second stage subproblem only has inequality constraints
        for (int scenario_index = 0; scenario_index < bi_dataset_size; ++scenario_index) {
            double tempCost = twoStageLP_secondStageCost_inequality(x_candidate, model_parameters.d, model_parameters.Di, model_parameters.Ci,bRHSDB.bi_database[0][scenario_index]);
            secondStageTotalCost += tempCost / bi_dataset_size;
            double tempTotalCost = tempCost + firstStageCost;
            double _n_ = scenario_index + 1;
            variance_P = variance_P * (_n_ - 1) / _n_ + tempTotalCost * tempTotalCost / _n_;
        }
    }
    else { // second stage subproblem is unconstrained
        std::cout << "Warning: Second stage subproblem is unconstrained! Return 0." << std::endl;
    }
    double _n_ = max(bi_dataset_size,be_dataset_size);
    validationResult result;
    result.mean = firstStageCost + secondStageTotalCost;
    result.variance = (variance_P - result.mean * result.mean) * _n_ / (_n_ - 1);
    result.num_dataPoint = max(bi_dataset_size,be_dataset_size);
    double halfMargin = result.Zalpha * sqrt(result.variance / _n_);
    result.CI_lower = result.mean - halfMargin;
    result.CI_upper = result.mean + halfMargin;
    // write computational results
    std::string outputResults_path = folder_path + "/validationResultsV2.txt";
    const char* outputResults_path_const = outputResults_path.c_str();
    std::fstream writeFile;
    writeFile.open(outputResults_path_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "***************************************************\n";
    writeFile << "Estimating the quality of candidate solution\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution             : ";
    long x_size = x_candidate.size();
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_candidate[index] << ", ";
    }
    writeFile << x_candidate[x_size - 1] << "\n";
    writeFile << "Number of data points          : " << result.num_dataPoint << "\n";
    writeFile << "Average validation cost        : " << result.mean << "\n";
    writeFile << "Variance                       : " << result.variance << "\n";
    writeFile << "Variance in estimating the mean: " << sqrt(result.variance/ _n_) << "\n";
    writeFile << result.alpha << "% confidence interval of expected cost: [" << result.CI_lower << ", " << result.CI_upper << "]\n";
    writeFile << "***************************************************\n";
    return result;
}

// find Dx
double Dx_estimateRectangle(const std::vector<double>&x, double lowerBound, double upperBound) {
    double Dx = 0;
    // size of x
    long x_size = x.size();
    // calculate Dx squared
    for (int x_index = 0; x_index < x_size; ++x_index) {
        Dx += max((x[x_index] - lowerBound) * (x[x_index] - lowerBound), (upperBound - x[x_index]) * (upperBound - x[x_index]));
    }
    Dx = sqrt(Dx);
    return Dx;
}

// test functions for debugging
void testTwoStageLPv1() {
    std::cout << "Test on Two Stage Stochastic Linear Programming Solver (NSQG)" << std::endl;
    std::string folder_path = "/Users/sonny/Documents/DebuggingPlayground/twoStageShipment02_EpanechnikovKernel2";
    std::vector<double> x(4,0);
    twoStageLP_random_b(folder_path, 8, 0.1, x);
}

void testTwoStageLPValidation() {
    std::cout << "Test on Obtaining Validation Cost in the Two Stage Stochastic Linear Programming" << std::endl;
    std::string folder_path = "/Users/sonny/Documents/DebuggingPlayground/twoStageShipment02_validation";
    std::vector<double> x(4,0);
    /*
     x[0] = 8.98;
     x[1] = 9.93;
     x[2] = 9.02;
     x[3] = 9.12;
     */
    double cost = twoStageLP_validation(folder_path, x);
    std::cout << "Candidate solution: (" << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << ")" << std::endl;
    std::cout << "Validation cost: " << cost << std::endl;
}

// test on projection
void projectionTest() {
    std::vector<double> x(4,-20);
    std::string folder_path = "/Users/sonny/Documents/numericalExperiment/twoStageShipment/case1/EpanechnikovKernel";
    std::string model_path = folder_path + "/model.txt";
    // read model
    twoStageParameters model_parameters = readTwoStageParameters(model_path);
    long A_rowsize = model_parameters.A.size();
    long A_colsize = model_parameters.A[0].size();
    std::vector<double> x_new = twoStageLP_projection(x, model_parameters.A, model_parameters.b, A_rowsize, A_colsize);
    std::cout << "Orignial point: (" << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << ")" << std::endl;
    std::cout << "Projected point: (" << x_new[0] << " " << x_new[1] << " " << x_new[2] << " " << x_new[3] << ")" << std::endl;
}

// iterative SAA
std::vector<double> twoStageLP_SAA(const std::string& folder_path) {
    // STEP 1: INITIALIZATION
    // create directory paths for database and model
    std::string be_DB_path = folder_path + "/be_DB.txt";
    std::string bi_DB_path = folder_path + "/bi_DB.txt";
    std::string model_path = folder_path + "/model.txt";
    std::string sto_path = folder_path + "/sto.txt";
    // convert all the paths into constant chars
    const char* be_DB_path_const = be_DB_path.c_str();
    const char* bi_DB_path_const = bi_DB_path.c_str();
    const char* model_path_const = model_path.c_str();
    const char* sto_path_const = sto_path.c_str();
    // create stream object
    std::ifstream readFile_be(be_DB_path_const);
    std::ifstream readFile_bi(bi_DB_path_const);
    // create database
    std::vector<std::vector<dataPoint>> be_DB;
    std::vector<std::vector<dataPoint>> bi_DB;
    // create model structure
    twoStageParameters model_parameters;
    // create sto object
    secondStageRHS bRHS;
    // read  be
    if (readFile_be.is_open()) {
        std::cout << "be_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read be database
        be_DB = readNonparametricDB(be_DB_path);
    }
    else {
        readFile_be.close(); // close the file
        std::cout << "be_DB is not found!" << std::endl;
    }
    // read bi
    if (readFile_bi.is_open()) {
        std::cout << "bi_DB is found." << std::endl;
        readFile_be.close(); // close the file
        // read bi database
        bi_DB = readNonparametricDB(bi_DB_path);
    }
    else {
        readFile_bi.close(); // close the file
        std::cout << "bi_DB is not found!" << std::endl;
    }
    // read model
    model_parameters = readTwoStageParameters(model_path);
    // read sto object
    bRHS = readStochastic(sto_path);
    // create be and bi in std::vector<double> form
    secondStageRHSDB bRHSDB = mergeDB_randomVector(be_DB, bi_DB, bRHS);
    // STEP 2: SOLVE TWO STAGE LINEAR PROGRAMMING
    long x_size = model_parameters.c.size();
    long y_size = model_parameters.d.size();
    long A_rowsize = model_parameters.A.size();
    long A_colsize = model_parameters.A[0].size();
    long be_numDatasets = bRHSDB.be_database.size();
    long bi_numDatasets = bRHSDB.bi_database.size();
    long numDatasets = 0;
    long numDatapoints = 0;
    if (be_numDatasets > 0) {
        numDatasets = be_numDatasets;
        numDatapoints = bRHSDB.be_database[0].size();
    }
    else {
        numDatasets = bi_numDatasets;
        numDatapoints = bRHSDB.bi_database[0].size();
    }
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray x(env,x_size,-IloInfinity,IloInfinity,ILOFLOAT); // first stage decision variable
    mod.add(x);
    std::vector<IloNumVarArray> y;// second stage decision variables
    for (int dataset_index = 0; dataset_index < numDatapoints; ++dataset_index) {
        IloNumVarArray y_temp(env,y_size,-IloInfinity,IloInfinity,ILOFLOAT); // second stage decision variable for one scenario
        y.push_back(y_temp);
        mod.add(y[dataset_index]);
    }
    // objective function
    IloExpr expr_obj(env);
    for (int x_index = 0; x_index < x_size; ++x_index) {
        if (model_parameters.c[x_index] != 0) {
            expr_obj += model_parameters.c[x_index] * x[x_index];
        }
    }
    double weight = 1.0 / ((double) numDatapoints);
    for (int dataPoint_index = 0; dataPoint_index < numDatapoints; ++dataPoint_index) {
        for (int y_index = 0; y_index < y_size; ++y_index) {
            if (model_parameters.d[y_index] != 0) {
                expr_obj += weight * model_parameters.d[y_index] * y[dataPoint_index][y_index];
            }
        }
    }
    IloObjective obj = IloMinimize(env, expr_obj);
    mod.add(obj);
    // constraints
    // first stage constraint
    for (int row_index = 0; row_index < A_rowsize; ++row_index) {
        IloExpr expr_cons1(env);
        for (int col_index = 0; col_index < A_colsize; ++col_index) {
            if (model_parameters.A[row_index][col_index] != 0) {
                expr_cons1 += model_parameters.A[row_index][col_index] * x[col_index];
            }
        }
        mod.add(expr_cons1 <= model_parameters.b[row_index]);
    }
    // second stage constraint
    // equality constraint
    if (be_numDatasets > 0) {
        for (int dataPoint_index = 0; dataPoint_index < numDatapoints; ++dataPoint_index) {
            long De_rowsize = model_parameters.De.size();
            for (int row_index = 0; row_index < De_rowsize; ++row_index) {
                IloExpr expr_consEq(env);
                for (int col_index = 0; col_index < y_size; ++col_index) {
                    if (model_parameters.De[row_index][col_index] != 0) {
                        expr_consEq += model_parameters.De[row_index][col_index] * y[dataPoint_index][col_index];
                    }
                }
                for (int col_index = 0; col_index < x_size; ++col_index) {
                    if (model_parameters.Ce[row_index][col_index] != 0) {
                        expr_consEq += model_parameters.Ce[row_index][col_index] * x[col_index];
                    }
                }
                mod.add(expr_consEq == bRHSDB.be_database[0][dataPoint_index][row_index]);
            }
        }
    }
    // inequality constraint
    if (bi_numDatasets > 0) {
        for (int dataPoint_index = 0; dataPoint_index < numDatapoints; ++dataPoint_index) {
            long Di_rowsize = model_parameters.Di.size();
            for (int row_index = 0; row_index < Di_rowsize; ++row_index) {
                IloExpr expr_consIneq(env);
                for (int col_index = 0; col_index < y_size; ++col_index) {
                    if (model_parameters.Di[row_index][col_index] != 0) {
                        expr_consIneq += model_parameters.Di[row_index][col_index] * y[dataPoint_index][col_index];
                    }
                }
                for (int col_index = 0; col_index < x_size; ++col_index) {
                    if (model_parameters.Ci[row_index][col_index] != 0) {
                        expr_consIneq += model_parameters.Ci[row_index][col_index] * x[col_index];
                    }
                }
                mod.add(expr_consIneq <= bRHSDB.bi_database[0][dataPoint_index][row_index]);
            }
        }
    }
    // solve the problem
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.solve();
    // get the estimated solution
    std::vector<double> x_est(x_size,0.0);
    for (int x_index = 0; x_index < x_size; ++x_index) {
        x_est[x_index] = cplex.getValue(x[x_index]);
    }
    // output the result
    std::string validation_outputPath = folder_path + "/validationResults.txt";
    const char* validation_outputPath_const = validation_outputPath.c_str();
    std::fstream writeFile;
    writeFile.open(validation_outputPath_const,std::fstream::app);
    // current time
    std::time_t currTime = std::time(nullptr);
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    std::cout << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Candidate solution   : ";
    std::cout << "Candidate solution   : ";
    for (int index = 0; index < x_size - 1; ++index) {
        writeFile << x_est[index] << ", ";
        std::cout << x_est[index] << ", ";
    }
    writeFile << x_est[x_size - 1] << "\n";
    std::cout << x_est[x_size - 1] << "\n";
    writeFile << "Optimal Cost         : " << cplex.getObjValue() << std::endl;
    std::cout << "Optimal Cost         : " << cplex.getObjValue() << std::endl;
    writeFile << "**********************************************************************\n";
    std::cout << "**********************************************************************\n";
    env.end();
    return x_est;
}
