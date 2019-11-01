//
//  nonparametricSQG.cpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/23/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#include "nonparametricSQG.hpp"
// test
void CplexConnectionTest(){ // test on whether it is successful to connect cplex
    std::cout << "Test on Cplex Connection" << std::endl;
    // solving min x^2
    IloEnv env;
    IloModel mod(env);
    IloNumVar x(env,-1,1,ILOFLOAT);
    mod.add(x);
    IloExpr expr(env);
    expr = x * x;
    IloObjective obj = IloMinimize(env,expr);
    mod.add(obj);
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setParam(IloCplex::RootAlg, IloCplex::Dual); // set up optimizer
    IloBool flag_solvable = cplex.solve();
    std::cout << "Final Results" << std::endl;
    std::cout << "Optimal solution is x = " << cplex.getValue(x) << std::endl;
    std::cout << "Optimal value is " << cplex.getObjValue() << std::endl;
    env.end();
}

void CplexExtremeRayTest() {
    std::cout << "Test on getting extreme rays" << std::endl;
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray x(env,2,0,IloInfinity,ILOFLOAT);
    mod.add(x);
    IloExpr expr1(env);
    expr1 = x[0] + x[1];
    IloObjective obj = IloMinimize(env,expr1);
    mod.add(obj);
    //IloConstraintArray constraints(env);
    IloRangeArray constraints(env);
    IloExpr expr2(env);
    expr2 = x[0] + x[1];
    constraints.add(expr2 <= -1);
    //IloRangeArray constraints2(env);
    IloExpr expr3(env);
    expr3 = x[0] - 0.5 * x[1];
    constraints.add(expr3 <= 0);
    IloExpr expr4(env);
    expr4 = x[0] - x[1];
    constraints.add(expr4 == 1);
    mod.add(constraints);
    //mod.add(constraints2);
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setParam(IloCplex::PreInd, false); // need to turn of presolve in order to get dual extreme rays
    cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);
    IloBool feasible_flag = cplex.solve();
    if (feasible_flag == IloFalse) {
        std::cout << "Problem is infeasible!" << std::endl;
    }
    IloNumArray y(env);
    //IloConstraintArray farkasConstraints(env);
    //cplex.getDuals(y, constraints);
    cplex.dualFarkas(constraints, y);
    //cplex.dualFarkas(farkasConstraints,y);
    std::cout << "y" << std::endl;
    std::cout << y[0] << std::endl;
    std::cout << y[1] << std::endl;
    std::cout << y[2] << std::endl;
    std::cout << y.getSize() << std::endl;
    double test_y = y[0];
    std::cout << test_y << std::endl;
    //IloNumArray z(env);
    //std::cout << "z" << std::endl;
    //cplex.dualFarkas(constraints2, z);
    //std::cout << z[1] << std::endl;
    env.end();
}

// test on getting extreme rays equality constraint
void CplexExtremeRayTest02() {
    std::cout << "Test on getting extreme rays" << std::endl;
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray x(env,2,0,IloInfinity,ILOFLOAT);
    mod.add(x);
    IloExpr expr1(env);
    expr1 = x[0] + x[1];
    IloObjective obj = IloMinimize(env,expr1);
    mod.add(obj);
    //IloConstraintArray constraints(env);
    IloRangeArray constraints(env);
    IloExpr expr2(env);
    expr2 = x[0] + x[1];
    constraints.add(expr2 == -1);
    //IloRangeArray constraints2(env);
    IloExpr expr3(env);
    expr3 = x[0] - x[1];
    constraints.add(expr3 <= 0);
    mod.add(constraints);
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setParam(IloCplex::PreInd, false); // need to turn of presolve in order to get dual extreme rays
    cplex.setParam(IloCplex::RootAlg, IloCplex::Dual);
    IloBool feasible_flag = cplex.solve();
    if (feasible_flag == IloFalse) {
        std::cout << "Problem is infeasible!" << std::endl;
    }
    IloNumArray y(env);
    //IloConstraintArray farkasConstraints(env);
    //cplex.getDuals(y, constraints);
    cplex.dualFarkas(constraints, y);
    //cplex.dualFarkas(farkasConstraints,y);
    std::cout << "y" << std::endl;
    std::cout << y[0] << std::endl;
    std::cout << y[1] << std::endl;
    std::cout << y.getSize() << std::endl;
    std::vector<double> pi;
    pi.push_back(-y[0]);
    pi.push_back(-y[1]);
    std::cout << "pi" << std::endl;
    std::cout << pi[0] << std::endl;
    std::cout << pi[1] << std::endl;
    env.end();
}

// two stage linear programming where only be and bi are random
// obtain dual multiplers of the second stage, given x (first stage decision variable), d, De, Ce, be, Di, Ci, bi
dualMultipliers twoStageLP_secondStageDual(const std::vector<double>& x, const std::vector<double>& d, const std::vector<std::vector<double>>& De, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be, const std::vector<std::vector<double>>& Di, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi){
    dualMultipliers pi; // initialize dual multipliers
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
    if (solvable_flag == IloTrue) {
        pi.feasible_flag = true; // tell the subproblem is feasible for given x, first stage decision variable
        if (equality_size > 0) {
            cplex.getDuals(dual_equality,constraintsEquality);
            for (int equality_index = 0; equality_index < equality_size; ++equality_index) {
                double pi_temp = -dual_equality[equality_index]; // sign of dual in cplex is opposite
                pi.equality.push_back(pi_temp);
            }
        }
        if (inequality_size > 0) {
            cplex.getDuals(dual_inequality,constraintsInequality);
            for (int inequality_index = 0; inequality_index < inequality_size; ++inequality_index) {
                double pi_temp = -dual_inequality[inequality_index];
                pi.inequality.push_back(pi_temp);
            }
        }
    }
    else {
        pi.feasible_flag = false; // tell the subproblem is infeasible for given x
    }
    env.end(); // end the environment
    return pi; // return the dual multipliers
}

dualMultipliers twoStageLP_secondStageDual_inequality(const std::vector<double>& x, const std::vector<double>& d,const std::vector<std::vector<double>>& Di, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi) {
    dualMultipliers pi; // initialize dual multipliers
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
    cplex.exportModel("/Users/sonny/Documents/DebuggingPlayground/twoStageShipment01/tss_subgradient_second_stage.lp");
    cplex.setOut(env.getNullStream());
    IloBool solvable_flag = cplex.solve();
    IloNumArray dual_equality(env);
    IloNumArray dual_inequality(env);
    if (solvable_flag == IloTrue) {
        pi.feasible_flag = true; // tell the subproblem is feasible for given x, first stage decision variable
        if (inequality_size > 0) {
            cplex.getDuals(dual_inequality,constraintsInequality);
            for (int inequality_index = 0; inequality_index < inequality_size; ++inequality_index) {
                double pi_temp = -dual_inequality[inequality_index];
                pi.inequality.push_back(pi_temp);
            }
        }
    }
    else {
        pi.feasible_flag = false; // tell the subproblem is infeasible for given x
    }
    env.end(); // end the environment
    return pi; // return the dual multipliers
} // end twoStageLP_secondStageDual without equality constraints

dualMultipliers twoStageLP_secondStageDual_equality(const std::vector<double>& x, const std::vector<double>& d,const std::vector<std::vector<double>>& De, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be){ // no inequality constraint
    dualMultipliers pi; // initialize dual multipliers
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
    if (solvable_flag == IloTrue) {
        pi.feasible_flag = true; // tell the subproblem is feasible for given x, first stage decision variable
        if (equality_size > 0) {
            cplex.getDuals(dual_equality,constraintsEquality);
            for (int equality_index = 0; equality_index < equality_size; ++equality_index) {
                double pi_temp = -dual_equality[equality_index]; // sign of dual in cplex is opposite
                pi.equality.push_back(pi_temp);
            }
        }
    }
    else {
        pi.feasible_flag = false; // tell the subproblem is infeasible for given x
    }
    env.end(); // end the environment
    return pi; // return the dual multipliers
} // end twoStageLP_secondStageDual without inequality constraints

// Nonparametric SQG solver for two stage linear programming
std::vector<double> twoStageLP_solver(twoStageParameters parameters,const std::vector<std::vector<std::vector<double>>>& be_database, const std::vector<std::vector<std::vector<double>>>& bi_database, const std::vector<std::vector<double>>& weight_database, double initial_stepsize, int max_iterates, std::vector<double> x_old) {
    // size of parameters
    long c_size = parameters.c.size();
    long A_rowsize = parameters.A.size();
    long A_colsize = 0;
    long De_rowsize = parameters.De.size();
    long Di_rowsize = parameters.Di.size();
    if (A_rowsize > 0) {
        A_colsize = parameters.A[0].size();
    }
    // initial variable
    std::vector<double> x_new(c_size,0.0);
    x_old = twoStageLP_projection(x_old, parameters.A, parameters.b, A_rowsize, A_colsize);
    // main loop
    for (int iterate_index = 0; iterate_index < max_iterates; ++iterate_index) {
        feasibilityCut feasibilityCut_scenario;
        bool feasibility_flag = true;
        std::vector<double> pi_equalityAGG(De_rowsize,0.0); // aggregate sum of duals
        std::vector<double> pi_inequalityAGG(Di_rowsize,0.0);
        // size of dataset (number of datapoints in current dataset)
        long be_dataset_size = 0;
        if ( iterate_index < be_database.size()) {
            be_dataset_size = be_database[iterate_index].size();
        }
        long bi_dataset_size = 0;
        if (iterate_index < bi_database.size()) {
            bi_dataset_size = bi_database[iterate_index].size();
        }
        long dataset_size = max(be_dataset_size,bi_dataset_size); // might exist some bugs
        // loop through all scenarios in this dataset
        for (int scenario_index = 0; scenario_index < dataset_size; ++scenario_index) {
            std::cout << "scenario_index: " << scenario_index << std::endl;
            dualMultipliers dualsTemp;
            if (be_dataset_size < 1 && bi_dataset_size > 0) { // no equality constraints in the second stage
                std::cout << "Second stage subproblem only has inequality constraints." << std::endl;
                dualsTemp = twoStageLP_secondStageDual_inequality(x_old, parameters.d, parameters.Di, parameters.Ci, bi_database[iterate_index][scenario_index]);
            }
            else if (be_dataset_size > 0 && bi_dataset_size < 1) { // no inequality constraints in the second stage
                std::cout << "Second stage subproblem only has equality constraints." << std::endl;
                dualsTemp = twoStageLP_secondStageDual_equality(x_old, parameters.d, parameters.De, parameters.Ce, be_database[iterate_index][scenario_index]);
            }
            else if (be_dataset_size > 0 && bi_dataset_size > 0) {
                std::cout << "Second stage subproblem has both equality constraints and inequality constraints." << std::endl;
                dualsTemp = twoStageLP_secondStageDual(x_old, parameters.d, parameters.De, parameters.Ce, be_database[iterate_index][scenario_index], parameters.Di, parameters.Ci, bi_database[iterate_index][scenario_index]); // obtain the dual multipliers of the second stage in current scenario
            }
            else { // second-stage feasible region is unconstrained
                std::cout << "Warning: Second-stage feasible region is unconstrained!" << std::endl;
                std::cout << "Warning: Return orignial estimate!" << std::endl;
                return x_old;
            }
            if (dualsTemp.feasible_flag == false) { // if second stage subproblem is infeasible
                // need to generate feasibility cut
                // obtain extrem ray of the subproblem
                dualMultipliers extremeRay_scenario;
                feasibilityCut feasibilityCut_scenario;
                if (be_dataset_size < 1 && bi_dataset_size > 0) { // no equality constraint in the subproblem
                    extremeRay_scenario = twoStageLP_secondStageExtremRay_inequality(x_old,  parameters.d, parameters.Di, parameters.Ci, bi_database[iterate_index][scenario_index]);
                    feasibilityCut_scenario = twoStageLP_feasibilityCutGeneration_inequality(extremeRay_scenario, parameters.Ci, bi_database[iterate_index][scenario_index]);
                }
                else if (be_dataset_size > 0 && bi_dataset_size < 1) { // no inequality constraint in the subproblem
                    extremeRay_scenario = twoStageLP_secondStageExtremRay_equality(x_old,  parameters.d, parameters.De, parameters.Ce, be_database[iterate_index][scenario_index]);
                    feasibilityCut_scenario = twoStageLP_feasibilityCutGeneration_equality(extremeRay_scenario, parameters.Ce, be_database[iterate_index][scenario_index]);
                }
                else if (be_dataset_size > 0 && bi_dataset_size > 0) {
                    extremeRay_scenario = twoStageLP_secondStageExtremRay(x_old,  parameters.d, parameters.De, parameters.Ce, be_database[iterate_index][scenario_index], parameters.Di, parameters.Ci, bi_database[iterate_index][scenario_index]);
                    feasibilityCut_scenario = twoStageLP_feasibilityCutGeneration(extremeRay_scenario, parameters.Ce, be_database[iterate_index][scenario_index], parameters.Ci, bi_database[iterate_index][scenario_index]);
                }
                else { // second-stage feasible region is unconstrained
                    std::cout << "Warning: Second-stage feasible region is unconstrained!" << std::endl;
                    std::cout << "Warning: Return orignial estimate!" << std::endl;
                    return x_old;
                }
                feasibility_flag = false; // problem is infeasible, set the flag to false
                break; // break the loop
            }
            for (int De_index = 0; De_index < De_rowsize; ++De_index) {
                pi_equalityAGG[De_index] += dualsTemp.equality[De_index] * weight_database[iterate_index][scenario_index];
            } // end for
            for (int Di_index = 0; Di_index < Di_rowsize; ++Di_index) {
                pi_inequalityAGG[Di_index] += dualsTemp.inequality[Di_index] * weight_database[iterate_index][scenario_index];
            }
        }// end for loop through all scenarios in the dataset
        // aggregated subgradient calculation
        std::vector<double> CeT_times_pi_equality(c_size,0.0);
        std::vector<double> CiT_times_pi_inequality(c_size,0.0);
        std::vector<double> Gx(c_size,0.0); // subgradient
        // estimated solution update
        double stepsize = initial_stepsize / ((double)(1 + iterate_index));
        if (feasibility_flag == true) {
            std::cout << "All the subproblems are feasible in the current iteration" << std::endl;
            for (int x_index = 0; x_index < c_size; ++x_index) {
                for (int pi_equality_index = 0; pi_equality_index < De_rowsize; ++pi_equality_index) {
                    CeT_times_pi_equality[x_index] += parameters.Ce[pi_equality_index][x_index] * pi_equalityAGG[pi_equality_index];
                }
                for (int pi_inequality_index = 0; pi_inequality_index < Di_rowsize; ++pi_inequality_index) {
                    CiT_times_pi_inequality[x_index] += parameters.Ci[pi_inequality_index][x_index] * pi_inequalityAGG[pi_inequality_index];
                }
                Gx[x_index] = CeT_times_pi_equality[x_index] + CiT_times_pi_inequality[x_index] + parameters.c[x_index];
            } // end for
            for (int x_index = 0; x_index < c_size; ++x_index) {
                x_old[x_index] = x_old[x_index] - Gx[x_index] * stepsize;
            }
        }
        else { // add new feasibility cut
            parameters.A.push_back(feasibilityCut_scenario.A_newRow);
            parameters.b.push_back(feasibilityCut_scenario.b_newRow);
            A_rowsize++; // add one row
        }
        // projection
        x_new = twoStageLP_projection(x_old, parameters.A, parameters.b, A_rowsize, A_colsize);
        std::cout << "Iteration: " << iterate_index << std::endl;
        std::cout << "x: (";
        // update x_old and go to the next iteration
        for (int x_index = 0; x_index < c_size; ++x_index) {
            x_old[x_index] = x_new[x_index];
            std::cout << x_new[x_index] << " ";
        }
        std::cout << ")" <<std::endl;
        std::cout << "**************************" << std::endl;
    }
    return  x_new;
}

// Nonparametric SQG solver which has output ability
std::vector<double> twoStageLP_solver(twoStageParameters parameters,const std::vector<std::vector<std::vector<double>>>& be_database, const std::vector<std::vector<std::vector<double>>>& bi_database, const std::vector<std::vector<double>>& weight_database, double initial_stepsize, int max_iterates, std::vector<double> x_old, const std::string& resultsOutput_path) {
    // timer
    std::clock_t time_start;
    time_start = std::clock();
    // current time
    std::time_t currTime = std::time(nullptr);
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // write initial setup
    writeFile << "*******************************************\n";
    writeFile << "Nonparametric Stochastic Quasi-Gradient Method\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Initial solution: ";
    long x_old_size = x_old.size();
    for (int index = 0; index < x_old_size - 1; ++index) {
        writeFile << x_old[index] << ", ";
    }
    writeFile << x_old[x_old_size - 1] << "\n";
    writeFile << "Initial stepsize: " << initial_stepsize << "\n";
    writeFile << "Max number of iterations: " << max_iterates << "\n";
    // size of parameters
    long c_size = parameters.c.size();
    long A_rowsize = parameters.A.size();
    long A_colsize = 0;
    long De_rowsize = parameters.De.size();
    long Di_rowsize = parameters.Di.size();
    if (A_rowsize > 0) {
        A_colsize = parameters.A[0].size();
    }
    // initial variable
    std::vector<double> x_new(c_size,0.0);
    x_old = twoStageLP_projection(x_old, parameters.A, parameters.b, A_rowsize, A_colsize);
    for (int x_index = 0; x_index < x_old_size; ++x_index) {
        x_new[x_index] = x_old[x_index];
    }
    // write main loop
    writeFile << "Main loop\n";
    // main loop
    for (int iterate_index = 0; iterate_index < max_iterates; ++iterate_index) {
        feasibilityCut feasibilityCut_scenario;
        bool feasibility_flag = true;
        std::vector<double> pi_equalityAGG(De_rowsize,0.0); // aggregate sum of duals
        std::vector<double> pi_inequalityAGG(Di_rowsize,0.0);
        // size of dataset (number of datapoints in current dataset)
        long be_dataset_size = 0;
        if ( iterate_index < be_database.size()) {
            be_dataset_size = be_database[iterate_index].size();
        }
        long bi_dataset_size = 0;
        if (iterate_index < bi_database.size()) {
            bi_dataset_size = bi_database[iterate_index].size();
        }
        long dataset_size = max(be_dataset_size,bi_dataset_size); // might exist some bugs
        // write number of scenarios
        writeFile << "------------------------------------------------\n";
        writeFile << "Iteration: " << iterate_index + 1 << "\n";
        writeFile << "Number of scenarios: " << dataset_size << "\n";
        // loop through all scenarios in this dataset
        for (int scenario_index = 0; scenario_index < dataset_size; ++scenario_index) {
            std::cout << "scenario_index: " << scenario_index << std::endl;
            dualMultipliers dualsTemp;
            if (be_dataset_size < 1 && bi_dataset_size > 0) { // no equality constraints in the second stage
                std::cout << "Second stage subproblem only has inequality constraints." << std::endl;
                dualsTemp = twoStageLP_secondStageDual_inequality(x_old, parameters.d, parameters.Di, parameters.Ci, bi_database[iterate_index][scenario_index]);
            }
            else if (be_dataset_size > 0 && bi_dataset_size < 1) { // no inequality constraints in the second stage
                std::cout << "Second stage subproblem only has equality constraints." << std::endl;
                dualsTemp = twoStageLP_secondStageDual_equality(x_old, parameters.d, parameters.De, parameters.Ce, be_database[iterate_index][scenario_index]);
            }
            else if (be_dataset_size > 0 && bi_dataset_size > 0) {
                std::cout << "Second stage subproblem has both equality constraints and inequality constraints." << std::endl;
                dualsTemp = twoStageLP_secondStageDual(x_old, parameters.d, parameters.De, parameters.Ce, be_database[iterate_index][scenario_index], parameters.Di, parameters.Ci, bi_database[iterate_index][scenario_index]); // obtain the dual multipliers of the second stage in current scenario
            }
            else { // second-stage feasible region is unconstrained
                // write second stage problem is unconstrained
                writeFile << "Warning: Second-Stage feasible region is unconstrained!\n";
                writeFile << "Warning: Return original estimate!\n";
                std::cout << "Warning: Second-stage feasible region is unconstrained!" << std::endl;
                std::cout << "Warning: Return orignial estimate!" << std::endl;
                return x_old;
            }
            if (dualsTemp.feasible_flag == false) { // if second stage subproblem is infeasible
                // need to generate feasibility cut
                // obtain extrem ray of the subproblem
                // write one second stage problem is infeasible
                writeFile << "Scenario " << scenario_index << " is infeasible, generate extreme ray\n";
                dualMultipliers extremeRay_scenario;
                feasibilityCut feasibilityCut_scenario;
                if (be_dataset_size < 1 && bi_dataset_size > 0) { // no equality constraint in the subproblem
                    extremeRay_scenario = twoStageLP_secondStageExtremRay_inequality(x_old,  parameters.d, parameters.Di, parameters.Ci, bi_database[iterate_index][scenario_index]);
                    feasibilityCut_scenario = twoStageLP_feasibilityCutGeneration_inequality(extremeRay_scenario, parameters.Ci, bi_database[iterate_index][scenario_index]);
                }
                else if (be_dataset_size > 0 && bi_dataset_size < 1) { // no inequality constraint in the subproblem
                    extremeRay_scenario = twoStageLP_secondStageExtremRay_equality(x_old,  parameters.d, parameters.De, parameters.Ce, be_database[iterate_index][scenario_index]);
                    feasibilityCut_scenario = twoStageLP_feasibilityCutGeneration_equality(extremeRay_scenario, parameters.Ce, be_database[iterate_index][scenario_index]);
                }
                else if (be_dataset_size > 0 && bi_dataset_size > 0) {
                    extremeRay_scenario = twoStageLP_secondStageExtremRay(x_old,  parameters.d, parameters.De, parameters.Ce, be_database[iterate_index][scenario_index], parameters.Di, parameters.Ci, bi_database[iterate_index][scenario_index]);
                    feasibilityCut_scenario = twoStageLP_feasibilityCutGeneration(extremeRay_scenario, parameters.Ce, be_database[iterate_index][scenario_index], parameters.Ci, bi_database[iterate_index][scenario_index]);
                }
                else { // second-stage feasible region is unconstrained
                    // write second stage problem is unconstrained
                    writeFile << "Warning: Second-Stage feasible region is unconstrained!\n";
                    writeFile << "Warning: Return original estimate!\n";
                    std::cout << "Warning: Second-stage feasible region is unconstrained!" << std::endl;
                    std::cout << "Warning: Return orignial estimate!" << std::endl;
                    return x_old;
                }
                feasibility_flag = false; // problem is infeasible, set the flag to false
                break; // break the loop
            }
            for (int De_index = 0; De_index < De_rowsize; ++De_index) {
                pi_equalityAGG[De_index] += dualsTemp.equality[De_index] * weight_database[iterate_index][scenario_index];
            } // end for
            for (int Di_index = 0; Di_index < Di_rowsize; ++Di_index) {
                pi_inequalityAGG[Di_index] += dualsTemp.inequality[Di_index] * weight_database[iterate_index][scenario_index];
            }
        }// end for loop through all scenarios in the dataset
        // aggregated subgradient calculation
        std::vector<double> CeT_times_pi_equality(c_size,0.0);
        std::vector<double> CiT_times_pi_inequality(c_size,0.0);
        std::vector<double> Gx(c_size,0.0); // subgradient
        // estimated solution update
        double stepsize = initial_stepsize / ((double)(1 + iterate_index));
        if (feasibility_flag == true) {
            std::cout << "All the subproblems are feasible in the current iteration" << std::endl;
            // write feasibility condition
            writeFile << "All the subproblems are feasible in the current iteration\n";
            for (int x_index = 0; x_index < c_size; ++x_index) {
                for (int pi_equality_index = 0; pi_equality_index < De_rowsize; ++pi_equality_index) {
                    CeT_times_pi_equality[x_index] += parameters.Ce[pi_equality_index][x_index] * pi_equalityAGG[pi_equality_index];
                }
                for (int pi_inequality_index = 0; pi_inequality_index < Di_rowsize; ++pi_inequality_index) {
                    CiT_times_pi_inequality[x_index] += parameters.Ci[pi_inequality_index][x_index] * pi_inequalityAGG[pi_inequality_index];
                }
                Gx[x_index] = CeT_times_pi_equality[x_index] + CiT_times_pi_inequality[x_index] + parameters.c[x_index];
            } // end for
            // write stepzie
            writeFile << "Stepsize: " << stepsize << "\n";
            // write quasigradient
            writeFile << "Quasigradient: ";
            for (int Gx_index = 0; Gx_index < c_size - 1; ++Gx_index) {
                writeFile << Gx[Gx_index] << ", ";
            }
            writeFile << Gx[c_size - 1] << "\n";
            // update estimate
            for (int x_index = 0; x_index < c_size; ++x_index) {
                x_old[x_index] = x_old[x_index] - Gx[x_index] * stepsize;
            }
        }
        else { // add new feasibility cut
            parameters.A.push_back(feasibilityCut_scenario.A_newRow);
            parameters.b.push_back(feasibilityCut_scenario.b_newRow);
            A_rowsize++; // add one row
        }
        // projection
        x_new = twoStageLP_projection(x_old, parameters.A, parameters.b, A_rowsize, A_colsize);
        std::cout << "Iteration: " << iterate_index + 1 << std::endl;
        std::cout << "x: (";
        // update x_old and go to the next iteration
        for (int x_index = 0; x_index < c_size; ++x_index) {
            x_old[x_index] = x_new[x_index];
            std::cout << x_new[x_index] << " ";
        }
        std::cout << ")" <<std::endl;
        std::cout << "**************************" << std::endl;
        // write new estimate
        writeFile << "New estimated solution: ";
        for (int x_index = 0; x_index < c_size - 1; ++x_index) {
            writeFile << x_new[x_index] << ", ";
        }
        writeFile << x_new[c_size - 1] << "\n";
    }
    // write time elapsed
    double duration = (std::clock() - time_start ) / (double) CLOCKS_PER_SEC;
    writeFile << "------------------------------------------------\n";
    writeFile << "Time elapsed(secs) : " << duration << "\n";
    writeFile << "*******************************************\n";
    writeFile.close(); // close writeFile
    return  x_new;
}

// Nonparametric SQG solver for one stage piecewise linear programming which has output ability
std::vector<double> oneStagePiecewiseLP_solver(oneStageParameters parameters, const oneStagePiecewiseLinearObjectiveDB& parametersDB_objective, double initial_stepsize, int max_iterates, std::vector<double> x_old, const std::string& resultsOutput_path) {
    // timer
    std::clock_t time_start;
    time_start = std::clock();
    // current time
    std::time_t currTime = std::time(nullptr);
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // write initial setup
    writeFile << "*******************************************\n";
    writeFile << "One Stage Piecewise Linear Programming\n";
    writeFile << "Nonparametric Stochastic Quasi-Gradient Method\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Initial solution: ";
    long x_old_size = x_old.size();
    for (int index = 0; index < x_old_size - 1; ++index) {
        writeFile << x_old[index] << ", ";
    }
    writeFile << x_old[x_old_size - 1] << "\n";
    writeFile << "Initial stepsize: " << initial_stepsize << "\n";
    writeFile << "Max number of iterations: " << max_iterates << "\n";
    // size of parameters
    long x_size = parameters.x.size();
    long A_rowsize = parameters.A.size();
    // initial solution and projection
    std::vector<double> x_new(x_size,0.0);
    x_old = LP_projection(x_old, parameters, A_rowsize, x_size); // projection
    // write main loop
    writeFile << "Main loop\n";
    // main loop
    for (int iterate_index = 0; iterate_index < max_iterates; ++iterate_index) {
        // number of scenarios in the current dataset
        long num_scenarios = parametersDB_objective.weight[iterate_index].size();
        // quasigradient of linear component in the objective
        std::vector<double> quasigradient_linearFunction(x_size,0.0);
        // quasigradient of max function
        std::vector<double> quasigradient_maxFunction(x_size,0.0);
        // write number of scenarios
        writeFile << "------------------------------------------------\n";
        writeFile << "Iteration: " << iterate_index + 1 << "\n";
        writeFile << "Number of scenarios: " << num_scenarios << "\n";
        // calculate subgradient for each scenario in the current dataset
        for (int scenario_index = 0; scenario_index < num_scenarios; ++scenario_index) {
            // max(piece1, piece2)
            double piece1 = 0;
            double piece2 = 0;
            // calculate piece1
            for (auto x_iterator : parametersDB_objective.c2[iterate_index][scenario_index]) {
                int x_index = x_iterator.first; // index of x
                double coefficient = x_iterator.second; // coeffcient for x[x_index]
                piece1 += coefficient * x_old[x_index];
            }
            piece1 += parametersDB_objective.e2[iterate_index][scenario_index];
            // calculate piece2
            for (auto x_iterator : parametersDB_objective.c3[iterate_index][scenario_index]) {
                int x_index = x_iterator.first; // index of x
                double coefficient = x_iterator.second; // coefficient for x[x_index]
                piece2 += coefficient * x_old[x_index];
            }
            piece2 += parametersDB_objective.e3[iterate_index][scenario_index];
            // compare the values of piece1 and piece2
            // and update quaisgradient of the weighted sum of max functions
            if (piece1 >= piece2) { // piece1 is active (larger)
                for (auto x_iterator : parametersDB_objective.c2[iterate_index][scenario_index]) {
                    int x_index = x_iterator.first; // index of x
                    double coefficient = x_iterator.second; // coeffcient for x[x_index]
                    quasigradient_maxFunction[x_index] += coefficient * parametersDB_objective.weight[iterate_index][scenario_index];
                } // end for
            } // piece2 is active (larger)
            else {
                for (auto x_iterator : parametersDB_objective.c3[iterate_index][scenario_index]) {
                    int x_index = x_iterator.first; // index of x
                    double coefficient = x_iterator.second; // coefficient for x[x_index]
                    quasigradient_maxFunction[x_index] += coefficient * parametersDB_objective.weight[iterate_index][scenario_index];
                }
            }
            for (auto x_iterator : parametersDB_objective.c1[iterate_index][scenario_index]) {
                int x_index = x_iterator.first; // index of x
                double coefficient = x_iterator.second; // coeffcient for x[x_index]
                quasigradient_linearFunction[x_index] += coefficient * parametersDB_objective.weight[iterate_index][scenario_index];
                //std::cout << "Index: " << x_index << std::endl;
                //std::cout << "coefficient : " << coefficient << std::endl;
                //std::cout << "Weight: " << parametersDB_objective.weight[iterate_index][scenario_index] << std::endl;
            }
        } // end for (each scenario)
        // step size
        double stepsize = initial_stepsize / ((double) iterate_index + 1.0);
        // update x_new
        for (int x_index = 0; x_index < x_size; ++x_index) {
            //std::cout << "Coefficient before piecewise linear function: " << parametersDB_objective.cPiecewise << std::endl;
            x_new[x_index] = x_old[x_index] - stepsize * (quasigradient_linearFunction[x_index] + parametersDB_objective.cPiecewise * quasigradient_maxFunction[x_index]);
            //std::cout << "x_old[" << x_index << "] = " << x_old[x_index] << std::endl;
            //std::cout << "x_new[" << x_index << "] = " << x_new[x_index] << std::endl;
        }
        // projection
        x_new = LP_projection(x_new, parameters, A_rowsize, x_size);
        // print out results
        std::cout << "Iteration: " << iterate_index + 1 << std::endl;
        std::cout << "x: (";
        // write stepzie
        writeFile << "Stepsize: " << stepsize << "\n";
        // write quasigradient
        writeFile << "Quasigradient: ";
        for (int x_index = 0; x_index < x_size - 1; ++x_index) {
            writeFile << quasigradient_linearFunction[x_index] + parametersDB_objective.cPiecewise * quasigradient_maxFunction[x_index] << ", ";
        }
        writeFile << quasigradient_linearFunction[x_size - 1] + parametersDB_objective.cPiecewise * quasigradient_maxFunction[x_size - 1] << "\n";
        // update x_old and write result
        writeFile << "New estimated solution: ";
        for (int x_index = 0; x_index < x_size - 1; ++x_index) {
            x_old[x_index] = x_new[x_index];
            std::cout << x_new[x_index] << ", ";
            writeFile << x_new[x_index] << ", ";
        }
        std::cout << x_new[x_size - 1] << ")\n";
        writeFile << x_new[x_size - 1] << "\n";
        std::cout << "**************************" << std::endl;
    }
    // write time elapsed
    double duration = (std::clock() - time_start ) / (double) CLOCKS_PER_SEC;
    writeFile << "------------------------------------------------\n";
    for (int x_index = 0; x_index < x_size - 1; ++x_index) {
        writeFile << "x[" << x_index << "] = " << x_new[x_index] << " Name: " << parameters.x[x_index] << "\n";
    }
    writeFile << "x[" << x_size - 1 << "] = " << x_new[x_size - 1] << " Name: " << parameters.x[x_size - 1] << "\n";
    writeFile << "Time elapsed(secs) : " << duration << "\n";
    writeFile << "*******************************************\n";
    writeFile.close(); // close writeFile
    return x_new;
} // end oneStagePiecewiseLP_solver

// Robust Nonparametric SQG solver which has output ability
std::vector<double> robustTwoStageLP_solver(twoStageParameters parameters,const std::vector<std::vector<std::vector<double>>>& be_database, const std::vector<std::vector<std::vector<double>>>& bi_database, const std::vector<std::vector<double>>& weight_database, int max_iterates, std::vector<double> x_old, double Dx, int m, double M, const std::string& resultsOutput_path) {
    // timer
    std::clock_t time_start;
    time_start = std::clock();
    // current time
    std::time_t currTime = std::time(nullptr);
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // write initial setup
    writeFile << "*******************************************\n";
    writeFile << "Robust Nonparametric Stochastic Quasi-Gradient Method\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Initial solution: ";
    long x_old_size = x_old.size();
    for (int index = 0; index < x_old_size - 1; ++index) {
        writeFile << x_old[index] << ", ";
    }
    writeFile << x_old[x_old_size - 1] << "\n";
    writeFile << "Max number of iterations: " << max_iterates << "\n";
    // size of parameters
    long c_size = parameters.c.size();
    long A_rowsize = parameters.A.size();
    long A_colsize = 0;
    long De_rowsize = parameters.De.size();
    long Di_rowsize = parameters.Di.size();
    if (A_rowsize > 0) {
        A_colsize = parameters.A[0].size();
    }
    // initial variable
    std::vector<double> x_new(c_size,0.0);
    std::vector<double> x_robust(c_size,0.0); // x_robust (take the average of the subsequence of estimates)
    x_old = twoStageLP_projection(x_old, parameters.A, parameters.b, A_rowsize, A_colsize);
    for (int x_index = 0; x_index < x_old_size; ++x_index) {
        x_new[x_index] = x_old[x_index];
        x_robust[x_index] = x_old[x_index];
    }
    // write main loop
    writeFile << "Main loop\n";
    // iterate index
    int iterate_index = 0;
    // main loop (outer loop)
    for (int outer_index = 0; outer_index < max_iterates; ++outer_index) {
        // initialization
        for (int x_index = 0; x_index < c_size; ++x_index) {
            x_robust[x_index] = 0;
        }
        // max number iterates in the inner loop
        int maxInner_iterates = (outer_index + 1) * m;
        std::cout << "Outer Loop: Iteration " << outer_index + 1 << std::endl;
        writeFile << "##################################################\n";
        writeFile << "Outer Loop: Iteration " << outer_index + 1 << "\n";
        // step size
        double stepsize = Dx / (M * sqrt((double) maxInner_iterates));
        for (int inner_index = 0; inner_index < maxInner_iterates; ++inner_index) {
            // compute quasigradient
            feasibilityCut feasibilityCut_scenario;
            bool feasibility_flag = true;
            std::vector<double> pi_equalityAGG(De_rowsize,0.0); // aggregate sum of duals
            std::vector<double> pi_inequalityAGG(Di_rowsize,0.0);
            // size of dataset (number of datapoints in current dataset)
            long be_dataset_size = 0;
            if ( iterate_index < be_database.size()) {
                be_dataset_size = be_database[iterate_index].size();
            }
            long bi_dataset_size = 0;
            if (iterate_index < bi_database.size()) {
                bi_dataset_size = bi_database[iterate_index].size();
            }
            long dataset_size = max(be_dataset_size,bi_dataset_size); // might exist some bugs
            // write number of scenarios
            writeFile << "------------------------------------------------\n";
            writeFile << "Inner loop: Iteration  " << inner_index + 1 << "\n";
            writeFile << "Number of scenarios: " << dataset_size << "\n";
            // loop through all scenarios in this dataset
            for (int scenario_index = 0; scenario_index < dataset_size; ++scenario_index) {
                std::cout << "scenario_index: " << scenario_index << std::endl;
                dualMultipliers dualsTemp;
                if (be_dataset_size < 1 && bi_dataset_size > 0) { // no equality constraints in the second stage
                    std::cout << "Second stage subproblem only has inequality constraints." << std::endl;
                    dualsTemp = twoStageLP_secondStageDual_inequality(x_old, parameters.d, parameters.Di, parameters.Ci, bi_database[iterate_index][scenario_index]);
                }
                else if (be_dataset_size > 0 && bi_dataset_size < 1) { // no inequality constraints in the second stage
                    std::cout << "Second stage subproblem only has equality constraints." << std::endl;
                    dualsTemp = twoStageLP_secondStageDual_equality(x_old, parameters.d, parameters.De, parameters.Ce, be_database[iterate_index][scenario_index]);
                }
                else if (be_dataset_size > 0 && bi_dataset_size > 0) {
                    std::cout << "Second stage subproblem has both equality constraints and inequality constraints." << std::endl;
                    dualsTemp = twoStageLP_secondStageDual(x_old, parameters.d, parameters.De, parameters.Ce, be_database[iterate_index][scenario_index], parameters.Di, parameters.Ci, bi_database[iterate_index][scenario_index]); // obtain the dual multipliers of the second stage in current scenario
                }
                else { // second-stage feasible region is unconstrained
                    // write second stage problem is unconstrained
                    writeFile << "Warning: Second-Stage feasible region is unconstrained!\n";
                    writeFile << "Warning: Return original estimate!\n";
                    std::cout << "Warning: Second-stage feasible region is unconstrained!" << std::endl;
                    std::cout << "Warning: Return orignial estimate!" << std::endl;
                    return x_old;
                }
                if (dualsTemp.feasible_flag == false) { // if second stage subproblem is infeasible
                    // need to generate feasibility cut
                    // obtain extrem ray of the subproblem
                    // write one second stage problem is infeasible
                    writeFile << "Scenario " << scenario_index << " is infeasible, generate extreme ray\n";
                    dualMultipliers extremeRay_scenario;
                    feasibilityCut feasibilityCut_scenario;
                    if (be_dataset_size < 1 && bi_dataset_size > 0) { // no equality constraint in the subproblem
                        extremeRay_scenario = twoStageLP_secondStageExtremRay_inequality(x_old,  parameters.d, parameters.Di, parameters.Ci, bi_database[iterate_index][scenario_index]);
                        feasibilityCut_scenario = twoStageLP_feasibilityCutGeneration_inequality(extremeRay_scenario, parameters.Ci, bi_database[iterate_index][scenario_index]);
                    }
                    else if (be_dataset_size > 0 && bi_dataset_size < 1) { // no inequality constraint in the subproblem
                        extremeRay_scenario = twoStageLP_secondStageExtremRay_equality(x_old,  parameters.d, parameters.De, parameters.Ce, be_database[iterate_index][scenario_index]);
                        feasibilityCut_scenario = twoStageLP_feasibilityCutGeneration_equality(extremeRay_scenario, parameters.Ce, be_database[iterate_index][scenario_index]);
                    }
                    else if (be_dataset_size > 0 && bi_dataset_size > 0) {
                        extremeRay_scenario = twoStageLP_secondStageExtremRay(x_old,  parameters.d, parameters.De, parameters.Ce, be_database[iterate_index][scenario_index], parameters.Di, parameters.Ci, bi_database[iterate_index][scenario_index]);
                        feasibilityCut_scenario = twoStageLP_feasibilityCutGeneration(extremeRay_scenario, parameters.Ce, be_database[iterate_index][scenario_index], parameters.Ci, bi_database[iterate_index][scenario_index]);
                    }
                    else { // second-stage feasible region is unconstrained
                        // write second stage problem is unconstrained
                        writeFile << "Warning: Second-Stage feasible region is unconstrained!\n";
                        writeFile << "Warning: Return original estimate!\n";
                        std::cout << "Warning: Second-stage feasible region is unconstrained!" << std::endl;
                        std::cout << "Warning: Return orignial estimate!" << std::endl;
                        return x_old;
                    }
                    feasibility_flag = false; // problem is infeasible, set the flag to false
                    break; // break the loop
                }
                for (int De_index = 0; De_index < De_rowsize; ++De_index) {
                    pi_equalityAGG[De_index] += dualsTemp.equality[De_index] * weight_database[iterate_index][scenario_index];
                } // end for
                for (int Di_index = 0; Di_index < Di_rowsize; ++Di_index) {
                    pi_inequalityAGG[Di_index] += dualsTemp.inequality[Di_index] * weight_database[iterate_index][scenario_index];
                }
            }// end for loop through all scenarios in the dataset
            // aggregated quasi-gradient calculation
            std::vector<double> CeT_times_pi_equality(c_size,0.0);
            std::vector<double> CiT_times_pi_inequality(c_size,0.0);
            std::vector<double> Gx(c_size,0.0); // subgradient
            // estimated solution update
            if (feasibility_flag == true) {
                std::cout << "All the subproblems are feasible in the current iteration" << std::endl;
                // write feasibility condition
                writeFile << "All the subproblems are feasible in the current iteration\n";
                for (int x_index = 0; x_index < c_size; ++x_index) {
                    for (int pi_equality_index = 0; pi_equality_index < De_rowsize; ++pi_equality_index) {
                        CeT_times_pi_equality[x_index] += parameters.Ce[pi_equality_index][x_index] * pi_equalityAGG[pi_equality_index];
                    }
                    for (int pi_inequality_index = 0; pi_inequality_index < Di_rowsize; ++pi_inequality_index) {
                        CiT_times_pi_inequality[x_index] += parameters.Ci[pi_inequality_index][x_index] * pi_inequalityAGG[pi_inequality_index];
                    }
                    Gx[x_index] = CeT_times_pi_equality[x_index] + CiT_times_pi_inequality[x_index] + parameters.c[x_index];
                } // end for
                // write stepzie
                writeFile << "Stepsize: " << stepsize << "\n";
                // write quasigradient
                writeFile << "Quasigradient: ";
                for (int Gx_index = 0; Gx_index < c_size - 1; ++Gx_index) {
                    writeFile << Gx[Gx_index] << ", ";
                }
                writeFile << Gx[c_size - 1] << "\n";
                // update estimate
                for (int x_index = 0; x_index < c_size; ++x_index) {
                    x_old[x_index] = x_old[x_index] - Gx[x_index] * stepsize;
                }
            }
            else { // add new feasibility cut
                parameters.A.push_back(feasibilityCut_scenario.A_newRow);
                parameters.b.push_back(feasibilityCut_scenario.b_newRow);
                A_rowsize++; // add one row
            }
            // projection
            x_new = twoStageLP_projection(x_old, parameters.A, parameters.b, A_rowsize, A_colsize);
            std::cout << "Iteration: " << inner_index + 1 << std::endl;
            std::cout << "x: (";
            // update x_old and go to the next iteration
            for (int x_index = 0; x_index < c_size; ++x_index) {
                x_old[x_index] = x_new[x_index];
                std::cout << x_new[x_index] << " ";
            }
            std::cout << ")" <<std::endl;
            std::cout << "**************************" << std::endl;
            writeFile << "New estimated solution in the subsequence: ";
            for (int x_index = 0; x_index < c_size - 1; ++x_index) {
                writeFile << x_new[x_index] << ", ";
                // calculate x_robust (take the average of the subsequence of estimates)
                x_robust[x_index] +=  1.0 / ((double) maxInner_iterates) * x_new[x_index];
            }
            writeFile << x_new[c_size - 1] << "\n";
            x_robust[c_size - 1] +=  1.0 / ((double) maxInner_iterates) * x_new[c_size - 1];
            writeFile << "------------------------------------------------\n";
            iterate_index++; // increment iterate
        }
        // write new estimate
        std::cout <<  "New estimated solution by Robust SQG: ";
        writeFile << "New estimated solution by Robust SQG: ";
        for (int x_index = 0; x_index < c_size - 1; ++x_index) {
            std::cout << x_robust[x_index] << ", ";
            writeFile << x_robust[x_index] << ", ";
        }
        std::cout << x_robust[c_size - 1] << "\n";
        writeFile << x_robust[c_size - 1] << "\n";
        writeFile << "##################################################\n";
    }
    // write time elapsed
    double duration = (std::clock() - time_start ) / (double) CLOCKS_PER_SEC;
    writeFile << "Time elapsed(secs) : " << duration << "\n";
    writeFile << "*******************************************\n";
    writeFile.close(); // close writeFile
    return  x_robust;
}

// Robust Nonparametric SQG solver for one stage piecewise linear programming which has output ability
std::vector<double> robustOneStagePiecewiseLP_solver(oneStageParameters parameters, const oneStagePiecewiseLinearObjectiveDB& parametersDB_objective, int max_iterates, std::vector<double> x_old, double Dx, int m, double M, const std::string& resultsOutput_path) {
    // timer
    std::clock_t time_start;
    time_start = std::clock();
    // current time
    std::time_t currTime = std::time(nullptr);
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // write initial setup
    writeFile << "*******************************************\n";
    writeFile << "One Stage Piecewise Linear Programming\n";
    writeFile << "Robust Nonparametric Stochastic Quasi-Gradient Method\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Initial solution: ";
    long x_old_size = x_old.size();
    for (int index = 0; index < x_old_size - 1; ++index) {
        writeFile << x_old[index] << ", ";
    }
    writeFile << x_old[x_old_size - 1] << "\n";
    writeFile << "Max number of iterations: " << max_iterates << "\n";
    // size of parameters
    long x_size = parameters.x.size();
    long A_rowsize = parameters.A.size();
    // initial solution and projection
    std::vector<double> x_new(x_size,0.0);
    std::vector<double> x_robust(x_size,0.0);
    x_old = LP_projection(x_old, parameters, A_rowsize, x_size); // projection
    // write main loop
    writeFile << "Main loop\n";
    // iterate index
    int iterate_index = 0;
    // main loop (outer loop)
    for (int outer_index = 0; outer_index < max_iterates; ++outer_index) {
        // refresh x_robust
        for (int x_index = 0; x_index < x_size; ++x_index) {
            x_robust[x_index] = 0;
        } // end for
        // max number iterates in the inner loop
        int maxInner_iterates = (outer_index + 1) * m;
        std::cout << "Outer Loop: Iteration " << outer_index + 1 << std::endl;
        writeFile << "##################################################\n";
        writeFile << "Outer Loop: Iteration " << outer_index + 1 << "\n";
        // step size
        double stepsize = Dx / (M * sqrt((double) maxInner_iterates));
        // inner loop
        for (int inner_index = 0; inner_index < maxInner_iterates; ++inner_index) {
            // number of scenarios in the current dataset
            long num_scenarios = parametersDB_objective.weight[iterate_index].size();
            // subgradient of linear component in the objective
            std::vector<double> quasigradient_linearFunction(x_size,0.0);
            // subgradient of max function
            std::vector<double> quasigradient_maxFunction(x_size,0.0);
            // write number of scenarios
            writeFile << "------------------------------------------------\n";
            writeFile << "Inner loop: Iteration  " << inner_index + 1 << "\n";
            writeFile << "Number of scenarios: " << num_scenarios << "\n";
            // calculate subgradient for each scenario in the current dataset
            for (int scenario_index = 0; scenario_index < num_scenarios; ++scenario_index) {
                // max(piece1, piece2)
                double piece1 = 0;
                double piece2 = 0;
                // calculate piece1
                for (auto x_iterator : parametersDB_objective.c2[iterate_index][scenario_index]) {
                    int x_index = x_iterator.first; // index of x
                    double coefficient = x_iterator.second; // coeffcient for x[x_index]
                    piece1 += coefficient * x_old[x_index];
                }
                piece1 += parametersDB_objective.e2[iterate_index][scenario_index];
                // calculate piece2
                for (auto x_iterator : parametersDB_objective.c3[iterate_index][scenario_index]) {
                    int x_index = x_iterator.first; // index of x
                    double coefficient = x_iterator.second; // coefficient for x[x_index]
                    piece2 += coefficient * x_old[x_index];
                }
                piece2 += parametersDB_objective.e3[iterate_index][scenario_index];
                // compare the values of piece1 and piece2
                // and update subgradient of the weighted sum of max functions
                if (piece1 >= piece2) { // piece1 is active
                    for (auto x_iterator : parametersDB_objective.c2[iterate_index][scenario_index]) {
                        int x_index = x_iterator.first; // index of x
                        double coefficient = x_iterator.second; // coeffcient for x[x_index]
                        quasigradient_maxFunction[x_index] += coefficient * parametersDB_objective.weight[iterate_index][scenario_index];
                    }// end for
                }
                else { // piece2 is larger
                    for (auto x_iterator : parametersDB_objective.c3[iterate_index][scenario_index]) {
                        int x_index = x_iterator.first; // index of x
                        double coefficient = x_iterator.second; // coefficient for x[x_index]
                        quasigradient_maxFunction[x_index] += coefficient * parametersDB_objective.weight[iterate_index][scenario_index];
                    }
                }
                for (auto x_iterator : parametersDB_objective.c1[iterate_index][scenario_index]) {
                    int x_index = x_iterator.first; // index of x
                    double coefficient = x_iterator.second; // coeffcient for x[x_index]
                    quasigradient_linearFunction[x_index] += coefficient * parametersDB_objective.weight[iterate_index][scenario_index];
                }
            } // end for (each scenario)
            // update x_new
            for (int x_index = 0; x_index < x_size; ++x_index) {
                //std::cout << "Coefficient before piecewise linear function: " << parametersDB_objective.cPiecewise << std::endl;
                x_new[x_index] = x_old[x_index] - stepsize * (quasigradient_linearFunction[x_index] + parametersDB_objective.cPiecewise * quasigradient_maxFunction[x_index]);
                //std::cout << "x_old[" << x_index << "] = " << x_old[x_index] << std::endl;
                //std::cout << "x_new[" << x_index << "] = " << x_new[x_index] << std::endl;
            }
            // projection
            x_new = LP_projection(x_new, parameters, A_rowsize, x_size);
            // write stepsize
            writeFile << "Stepsize: " << stepsize << "\n";
            // write quasigradient
            writeFile << "Quasigradient: ";
            for (int x_index = 0; x_index < x_size - 1; ++x_index) {
                writeFile << quasigradient_linearFunction[x_index] + parametersDB_objective.cPiecewise * quasigradient_maxFunction[x_index] << ", ";
            }
            writeFile << quasigradient_linearFunction[x_size - 1] + parametersDB_objective.cPiecewise * quasigradient_maxFunction[x_size - 1] << "\n";
            // update x_old and write result
            std::cout << "Iteration: " << inner_index + 1 << std::endl;
            std::cout << "x: (";
            writeFile << "New estimated solution in the subsequence: ";
            for (int x_index = 0; x_index < x_size - 1; ++x_index) {
                x_old[x_index] = x_new[x_index];
                // calculate x_robust (take the average of the subsequence of estimates)
                x_robust[x_index] +=  1.0 / ((double) maxInner_iterates) * x_new[x_index];
                std::cout << x_new[x_index] << ", ";
                writeFile << x_new[x_index] << ", ";
            }
            std::cout << x_new[x_size - 1] << ")\n";
            std::cout << "**************************" << std::endl;
            writeFile << x_new[x_size - 1] << "\n";
            x_old[x_size - 1] = x_new[x_size - 1]; // update x_old
            x_robust[x_size - 1] +=  1.0 / ((double) maxInner_iterates) * x_new[x_size - 1]; // update x_robust
            writeFile << "------------------------------------------------\n";
            iterate_index++; // increment iterate
        } // end for inner loop
        // write new estimate
        std::cout <<  "New estimated solution by Robust SQG: ";
        writeFile << "New estimated solution by Robust SQG: ";
        for (int x_index = 0; x_index < x_size - 1; ++x_index) {
            std::cout << x_robust[x_index] << ", ";
            writeFile << x_robust[x_index] << ", ";
        }
        std::cout << x_robust[x_size - 1] << "\n";
        writeFile << x_robust[x_size - 1] << "\n";
        writeFile << "##################################################\n";
        // under construction...
    } // end for main loop
    // write time elapsed
    double duration = (std::clock() - time_start ) / (double) CLOCKS_PER_SEC;
    writeFile << "------------------------------------------------\n";
    writeFile << "Estimated solution by Robust SQG: ";
    for (int x_index = 0; x_index < x_size - 1; ++x_index) {
        writeFile << "x[" << x_index << "] = " << x_robust[x_index] << " Name: " << parameters.x[x_index] << "\n";
    }
    writeFile << "x[" << x_size - 1 << "] = " << x_robust[x_size - 1] << " Name: " << parameters.x[x_size - 1] << "\n";
    writeFile << "Time elapsed(secs) : " << duration << "\n";
    writeFile << "*******************************************\n";
    writeFile.close(); // close writeFile
    return x_robust;
}

std::vector<double> robustOneStagePiecewiseLP_solverDEBUG(oneStageParameters parameters, const oneStagePiecewiseLinearObjectiveDB& parametersDB_objective, int max_iterates, std::vector<double> x_old, double Dx, int m, double M, const std::string& resultsOutput_path) {
    // timer
    std::clock_t time_start;
    time_start = std::clock();
    // current time
    std::time_t currTime = std::time(nullptr);
    // initialization of output file
    const char* writeFilePath = resultsOutput_path.c_str();
    std::fstream writeFile;
    writeFile.open(writeFilePath,std::fstream::app); // append results to the end of the file
    // write initial setup
    writeFile << "*******************************************\n";
    writeFile << "One Stage Piecewise Linear Programming\n";
    writeFile << "Robust Nonparametric Stochastic Quasi-Gradient Method\n";
    writeFile << std::put_time(localtime(&currTime), "%c %Z") << "\n"; // write current time
    writeFile << "Initial solution: ";
    long x_old_size = x_old.size();
    for (int index = 0; index < x_old_size - 1; ++index) {
        writeFile << x_old[index] << ", ";
    }
    writeFile << x_old[x_old_size - 1] << "\n";
    writeFile << "Max number of iterations: " << max_iterates << "\n";
    // size of parameters
    long x_size = parameters.x.size();
    long A_rowsize = parameters.A.size();
    // initial solution and projection
    std::vector<double> x_new(x_size,0.0);
    std::vector<double> x_robust(x_size,0.0);
    x_old = LP_projection(x_old, parameters, A_rowsize, x_size); // projection
    writeFile << "--DEBUG--\n";
    writeFile << "Projecting initial solution onto the feasible region" << std::endl;
    writeFile << "x: ";
    for (int x_index = 0; x_index < x_size - 1; ++x_index) {
        writeFile << x_old[x_index] << ", ";
    }
    writeFile << x_old[x_size - 1] << std::endl;
    writeFile << "--DEBUG--\n";
    // write main loop
    writeFile << "Main loop\n";
    // iterate index
    int iterate_index = 0;
    // main loop (outer loop)
    for (int outer_index = 0; outer_index < max_iterates; ++outer_index) {
        // refresh x_robust
        for (int x_index = 0; x_index < x_size; ++x_index) {
            x_robust[x_index] = 0;
        } // end for
        // max number iterates in the inner loop
        int maxInner_iterates = (outer_index + 1) * m;
        std::cout << "Outer Loop: Iteration " << outer_index + 1 << std::endl;
        writeFile << "##################################################\n";
        writeFile << "Outer Loop: Iteration " << outer_index + 1 << "\n";
        // step size
        double stepsize = Dx / (M * sqrt((double) maxInner_iterates));
        // inner loop
        for (int inner_index = 0; inner_index < maxInner_iterates; ++inner_index) {
            // number of scenarios in the current dataset
            long num_scenarios = parametersDB_objective.weight[iterate_index].size();
            // subgradient of linear component in the objective
            std::vector<double> quasigradient_linearFunction(x_size,0.0);
            // subgradient of max function
            std::vector<double> quasigradient_maxFunction(x_size,0.0);
            // write number of scenarios
            writeFile << "------------------------------------------------\n";
            writeFile << "Inner loop: Iteration  " << inner_index + 1 << "\n";
            writeFile << "Number of scenarios: " << num_scenarios << "\n";
            // calculate subgradient for each scenario in the current dataset
            for (int scenario_index = 0; scenario_index < num_scenarios; ++scenario_index) {
                // max(piece1, piece2)
                double piece1 = 0;
                double piece2 = 0;
                // calculate piece1
                for (auto x_iterator : parametersDB_objective.c2[iterate_index][scenario_index]) {
                    int x_index = x_iterator.first; // index of x
                    double coefficient = x_iterator.second; // coeffcient for x[x_index]
                    piece1 += coefficient * x_old[x_index];
                }
                piece1 += parametersDB_objective.e2[iterate_index][scenario_index];
                writeFile << "--DEBUG--\n";
                writeFile << "Piece 1: " << piece1 << std::endl;
                // calculate piece2
                for (auto x_iterator : parametersDB_objective.c3[iterate_index][scenario_index]) {
                    int x_index = x_iterator.first; // index of x
                    double coefficient = x_iterator.second; // coefficient for x[x_index]
                    piece2 += coefficient * x_old[x_index];
                }
                piece2 += parametersDB_objective.e3[iterate_index][scenario_index];
                writeFile << "Piece 2: " << piece2 << std::endl;
                writeFile << "--DEBUG--\n";
                // compare the values of piece1 and piece2
                // and update subgradient of the weighted sum of max functions
                if (piece1 >= piece2) { // piece1 is active
                    writeFile << "--DEBUG--\n";
                    writeFile << "Piece 1 is larger\n";
                    writeFile << "--DEBUG--\n";
                    for (auto x_iterator : parametersDB_objective.c2[iterate_index][scenario_index]) {
                        int x_index = x_iterator.first; // index of x
                        double coefficient = x_iterator.second; // coeffcient for x[x_index]
                        writeFile << "--DEBUG--\n";
                        writeFile << "index: " << x_index << " coefficient: " << coefficient << std::endl;
                        writeFile << "--DEBUG--\n";
                        quasigradient_maxFunction[x_index] += coefficient * parametersDB_objective.weight[iterate_index][scenario_index];
                    }// end for
                }
                else { // piece2 is larger
                    writeFile << "--DEBUG--\n";
                    writeFile << "Piece 2 is larger\n";
                    writeFile << "--DEBUG--\n";
                    for (auto x_iterator : parametersDB_objective.c3[iterate_index][scenario_index]) {
                        int x_index = x_iterator.first; // index of x
                        double coefficient = x_iterator.second; // coefficient for x[x_index]
                        writeFile << "--DEBUG--\n";
                        writeFile << "index: " << x_index << " coefficient: " << coefficient << std::endl;
                        writeFile << "--DEBUG--\n";
                        quasigradient_maxFunction[x_index] += coefficient * parametersDB_objective.weight[iterate_index][scenario_index];
                    }
                }
                for (auto x_iterator : parametersDB_objective.c1[iterate_index][scenario_index]) {
                    int x_index = x_iterator.first; // index of x
                    double coefficient = x_iterator.second; // coeffcient for x[x_index]
                    writeFile << "--DEBUG--\n";
                    writeFile << "c1\n";
                    writeFile << "index: " << x_index << " coefficient: " << coefficient << std::endl;
                    writeFile << "--DEBUG--\n";
                    quasigradient_linearFunction[x_index] += coefficient * parametersDB_objective.weight[iterate_index][scenario_index];
                }
            } // end for (each scenario)
            // update x_new
            for (int x_index = 0; x_index < x_size; ++x_index) {
                //std::cout << "Coefficient before piecewise linear function: " << parametersDB_objective.cPiecewise << std::endl;
                x_new[x_index] = x_old[x_index] - stepsize * (quasigradient_linearFunction[x_index] + parametersDB_objective.cPiecewise * quasigradient_maxFunction[x_index]);
                writeFile << "--DEBUG--\n";
                writeFile << "x_old[" << x_index << "] = " << x_old[x_index] << "\n";
                writeFile << "stepsize: " << stepsize << "\n";
                writeFile << "quasigradient_linearFunction[" << x_index << "] = " << quasigradient_linearFunction[x_index] << "\n";
                writeFile << "parametersDB_objective.cPiecewise: " << parametersDB_objective.cPiecewise << "\n";
                writeFile << "quasigradient_maxFunction[" << x_index << "] = " << quasigradient_maxFunction[x_index] << std::endl;
                writeFile << "x_new[" << x_index << "] = " << x_new[x_index] << "\n";
                writeFile << "--DEBUG--\n";
                //std::cout << "x_old[" << x_index << "] = " << x_old[x_index] << std::endl;
                //std::cout << "x_new[" << x_index << "] = " << x_new[x_index] << std::endl;
            }
            // projection
            x_new = LP_projection(x_new, parameters, A_rowsize, x_size);
            // write stepsize
            writeFile << "Stepsize: " << stepsize << "\n";
            // write quasigradient
            writeFile << "Quasigradient: ";
            for (int x_index = 0; x_index < x_size - 1; ++x_index) {
                writeFile << quasigradient_linearFunction[x_index] + parametersDB_objective.cPiecewise * quasigradient_maxFunction[x_index] << ", ";
            }
            writeFile << quasigradient_linearFunction[x_size - 1] + parametersDB_objective.cPiecewise * quasigradient_maxFunction[x_size - 1] << "\n";
            // update x_old and write result
            std::cout << "Iteration: " << inner_index + 1 << std::endl;
            std::cout << "x: (";
            writeFile << "New estimated solution in the subsequence: ";
            for (int x_index = 0; x_index < x_size - 1; ++x_index) {
                x_old[x_index] = x_new[x_index];
                // calculate x_robust (take the average of the subsequence of estimates)
                x_robust[x_index] +=  1.0 / ((double) maxInner_iterates) * x_new[x_index];
                std::cout << x_new[x_index] << ", ";
                writeFile << x_new[x_index] << ", ";
            }
            std::cout << x_new[x_size - 1] << ")\n";
            std::cout << "**************************" << std::endl;
            writeFile << x_new[x_size - 1] << "\n";
            x_old[x_size - 1] = x_new[x_size - 1]; // update x_old
            x_robust[x_size - 1] +=  1.0 / ((double) maxInner_iterates) * x_new[x_size - 1]; // update x_robust
            writeFile << "------------------------------------------------\n";
            iterate_index++; // increment iterate
        } // end for inner loop
        // write new estimate
        std::cout <<  "New estimated solution by Robust SQG: ";
        writeFile << "New estimated solution by Robust SQG: ";
        for (int x_index = 0; x_index < x_size - 1; ++x_index) {
            std::cout << x_robust[x_index] << ", ";
            writeFile << x_robust[x_index] << ", ";
        }
        std::cout << x_robust[x_size - 1] << "\n";
        writeFile << x_robust[x_size - 1] << "\n";
        writeFile << "##################################################\n";
        // under construction...
    } // end for main loop
    // write time elapsed
    double duration = (std::clock() - time_start ) / (double) CLOCKS_PER_SEC;
    writeFile << "------------------------------------------------\n";
    writeFile << "Estimated solution by Robust SQG: ";
    for (int x_index = 0; x_index < x_size - 1; ++x_index) {
        writeFile << "x[" << x_index << "] = " << x_robust[x_index] << " Name: " << parameters.x[x_index] << "\n";
    }
    writeFile << "x[" << x_size - 1 << "] = " << x_robust[x_size - 1] << " Name: " << parameters.x[x_size - 1] << "\n";
    writeFile << "Time elapsed(secs) : " << duration << "\n";
    writeFile << "*******************************************\n";
    writeFile.close(); // close writeFile
    return x_robust;
}

// twoStageLP projection
std::vector<double> twoStageLP_projection(const std::vector<double> x, const std::vector<std::vector<double>>& A, const std::vector<double> b, const long& A_rowsize, const long& A_colsize) {
    std::vector<double> x_projected(A_colsize,0.0);
    // solve a quadratic programming
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray x_temp(env,A_colsize,-IloInfinity,IloInfinity,ILOFLOAT);
    mod.add(x_temp);
    IloExpr expr_obj(env);
    for (int x_index = 0; x_index < A_colsize; ++x_index) {
        expr_obj += x_temp[x_index] * x_temp[x_index] - 2.0 * x_temp[x_index] * x[x_index];
    }
    IloObjective obj = IloMinimize(env,expr_obj); // objective function
    mod.add(obj);
    // constraints
    for (int A_rowIndex = 0; A_rowIndex < A_rowsize; ++A_rowIndex) {
        IloExpr expr_constraint(env);
        for (int x_index = 0; x_index < A_colsize; ++x_index) {
            expr_constraint += A[A_rowIndex][x_index] * x_temp[x_index];
        }
        mod.add(expr_constraint <= b[A_rowIndex]); // add the constraint to the model
    }
    // create cplex environment
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setOut(env.getNullStream());
    cplex.solve();
    // obtain the projected point
    for (int x_index = 0; x_index < A_colsize; ++x_index) {
        x_projected[x_index] = cplex.getValue(x_temp[x_index]);
        //std::cout << cplex.getValue(x_temp[x_index]) << std::endl;
    }
    env.end();
    return x_projected;
} // end twoStageLP projection

// projecting x onto a convex feasible region
std::vector<double> LP_projection(const std::vector<double>& x, const oneStageParameters& parameters, long A_rowsize, long x_size) {
    // solve a quadratic programming
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray x_temp(env,x_size,-IloInfinity,IloInfinity,ILOFLOAT);
    mod.add(x_temp);
    IloExpr expr_obj(env);
    // create objective
    for (int x_index = 0; x_index < x_size; ++x_index) {
        expr_obj += x_temp[x_index] * x_temp[x_index] - 2.0 * x_temp[x_index] * x[x_index];
    }
    IloObjective obj = IloMinimize(env, expr_obj);
    mod.add(obj);
    //std::cout << "size of x: " << x_size << std::endl;
    // constraints
    for (int row_index = 0; row_index < A_rowsize; ++row_index) {
        IloExpr expr_constrant(env);
        for (auto a_iterator : parameters.A[row_index]) {
            int x_index = a_iterator.first;
            double coefficient = a_iterator.second;
            //std::cout << "index of x: " << x_index << std::endl;
            expr_constrant += coefficient * x_temp[x_index];
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
    //cplex.exportModel("Users/sonny/Documents/S&P100/case1/kNN/clpexmodel");
    cplex.setOut(env.getNullStream()); // stop the output
    cplex.solve();
    // obtain the projected point
    std::vector<double> x_projected(x_size,0);
    for (int x_index = 0; x_index < x_size; ++x_index) {
        x_projected[x_index] = cplex.getValue(x_temp[x_index]);
    }
    env.end(); // end cplex environment
    return x_projected;
}

/* functions for feasibility cut generation
 */
dualMultipliers twoStageLP_secondStageExtremRay(const std::vector<double>& x, const std::vector<double>& d, const std::vector<std::vector<double>>& De, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be, const std::vector<std::vector<double>>& Di, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi) {
    std::cout << "Second stage subproblem is infeasible. Searching for extreme ray." << std::endl;
    dualMultipliers extremeRay;
    extremeRay.feasible_flag = false;
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
    long bi_rowsize = bi.size();
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
    IloRangeArray constraints(env);
    // equality constraints
    for (int equality_index = 0; equality_index < equality_size; ++equality_index) {
        IloExpr expr_equality(env);
        for (int De_index = 0; De_index < De_colsize; ++De_index) {
            expr_equality += De[equality_index][De_index] * y[De_index];
        }
        double Ce_times_x = Ce[equality_index] * x;
        expr_equality += Ce_times_x - be[equality_index];
        constraints.add(expr_equality == 0);
    }
    // inequality constraints
    for (int inequality_index = 0; inequality_index < inequality_size; ++inequality_index) {
        IloExpr expr_inequality(env);
        for (int Di_index = 0; Di_index < Di_colsize; ++Di_index) {
            expr_inequality += Di[inequality_index][Di_index] * y[Di_index];
        }
        double Ci_times_x = Ci[inequality_index] * x;
        expr_inequality += Ci_times_x - bi[inequality_index];
        constraints.add(expr_inequality <= 0);
    }
    mod.add(constraints); // add constraints to the model
    // set up cplex solver
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::PreInd, false); // need to turn of presolve in order to get dual extreme rays
    cplex.setParam(IloCplex::RootAlg, IloCplex::Dual); // use dual simplex optimizer
    IloBool solvable_flag = cplex.solve();
    IloNumArray twoStageLP_IloExtremeRay(env);
    cplex.dualFarkas(constraints, twoStageLP_IloExtremeRay);
    // first equality_size componenents are for equality constraints and the rest components are for inequality constraints
    long extremRay_size = twoStageLP_IloExtremeRay.getSize(); // get the size of extrem ray
    if (extremRay_size != equality_size + inequality_size) {
        std::cout << "Warning: Second Stage Problem Setup is wrong!" << std::endl;
    }
    for (int equality_index = 0; equality_index < equality_size; ++equality_index) {
        extremeRay.equality.push_back(-twoStageLP_IloExtremeRay[equality_index]);
    }
    for (int inequality_index = 0; inequality_index < inequality_size; ++inequality_index) {
        extremeRay.inequality.push_back(-twoStageLP_IloExtremeRay[equality_size + inequality_index]);
    }
    return extremeRay;
}

// there is no equality constraint in the second stage subproblem
dualMultipliers twoStageLP_secondStageExtremRay_inequality(const std::vector<double>& x, const std::vector<double>& d, const std::vector<std::vector<double>>& Di, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi) {
    std::cout << "Second stage subproblem is infeasible. Searching for extreme ray." << std::endl;
    dualMultipliers extremeRay;
    extremeRay.feasible_flag = false;
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
    long bi_rowsize = bi.size();
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
    IloRangeArray constraints(env);
    // inequality constraints
    for (int inequality_index = 0; inequality_index < inequality_size; ++inequality_index) {
        IloExpr expr_inequality(env);
        for (int Di_index = 0; Di_index < Di_colsize; ++Di_index) {
            expr_inequality += Di[inequality_index][Di_index] * y[Di_index];
        }
        double Ci_times_x = Ci[inequality_index] * x;
        expr_inequality += Ci_times_x - bi[inequality_index];
        constraints.add(expr_inequality <= 0);
    }
    mod.add(constraints); // add constraints to the model
    // set up cplex solver
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setParam(IloCplex::PreInd, false); // need to turn of presolve in order to get dual extreme rays
    cplex.setParam(IloCplex::RootAlg, IloCplex::Dual); // use dual simplex optimizer
    cplex.setOut(env.getNullStream());
    IloBool solvable_flag = cplex.solve();
    IloNumArray twoStageLP_IloExtremeRay(env);
    cplex.dualFarkas(constraints, twoStageLP_IloExtremeRay);
    // first equality_size componenents are for equality constraints and the rest components are for inequality constraints
    long extremRay_size = twoStageLP_IloExtremeRay.getSize(); // get the size of extrem ray
    if (extremRay_size != inequality_size) {
        std::cout << "Warning: Number of extreme rays does not match the number of inequality constraints!" << std::endl;
    }
    for (int inequality_index = 0; inequality_index < inequality_size; ++inequality_index) {
        extremeRay.inequality.push_back(-twoStageLP_IloExtremeRay[inequality_index]);
    }
    return extremeRay;
} // twoStageLP_secondStageExtremeRay_inequality

// there is no inequality constraint in the second stage subproblem
dualMultipliers twoStageLP_secondStageExtremRay_equality(const std::vector<double>& x, const std::vector<double>& d, const std::vector<std::vector<double>>& De, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be) {
    std::cout << "Second stage subproblem is infeasible. Searching for extreme ray." << std::endl;
    dualMultipliers extremeRay;
    extremeRay.feasible_flag = false;
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
    IloRangeArray constraints(env);
    // equality constraints
    for (int equality_index = 0; equality_index < equality_size; ++equality_index) {
        IloExpr expr_equality(env);
        for (int De_index = 0; De_index < De_colsize; ++De_index) {
            expr_equality += De[equality_index][De_index] * y[De_index];
        }
        double Ce_times_x = Ce[equality_index] * x;
        expr_equality += Ce_times_x - be[equality_index];
        constraints.add(expr_equality == 0);
    }
    mod.add(constraints); // add constraints to the model
    // set up cplex solver
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setParam(IloCplex::PreInd, false); // need to turn of presolve in order to get dual extreme rays
    cplex.setParam(IloCplex::RootAlg, IloCplex::Dual); // use dual simplex optimizer
    cplex.setOut(env.getNullStream());
    IloBool solvable_flag = cplex.solve();
    IloNumArray twoStageLP_IloExtremeRay(env);
    cplex.dualFarkas(constraints, twoStageLP_IloExtremeRay);
    // first equality_size componenents are for equality constraints and the rest components are for inequality constraints
    long extremRay_size = twoStageLP_IloExtremeRay.getSize(); // get the size of extrem ray
    if (extremRay_size != equality_size) {
        std::cout << "Warning: Number of extrem rays does not match the number of equality constraints!" << std::endl;
    }
    for (int equality_index = 0; equality_index < equality_size; ++equality_index) {
        extremeRay.equality.push_back(-twoStageLP_IloExtremeRay[equality_index]);
    }
    return extremeRay;
} // twoStageLP_secondStageExtremeRay_equality

// generate feasibility cut for given extreme ray
feasibilityCut twoStageLP_feasibilityCutGeneration(const dualMultipliers& extremeRay, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi) {
    feasibilityCut cut_scenario;
    // size of parameters
    long Ce_rowsize = Ce.size();
    long Ce_colsize = 0;
    if (Ce_rowsize > 0) {
        Ce_colsize = Ce[0].size();
    }
    long Ci_rowsize = Ci.size();
    long Ci_colsize = 0;
    if (Ci_rowsize > 0) {
        Ci_colsize = Ci[0].size();
    }
    long x_size = max(Ce_colsize, Ci_colsize); // size of x
    std::vector<double> a_new(x_size,0.0);
    for (int x_index = 0; x_index < x_size; ++x_index) {
        for (int equality_index = 0; equality_index < Ce_rowsize; ++equality_index) {
            a_new[x_index] += Ce[equality_index][x_index] * extremeRay.equality[equality_index];
        }
        for (int inequality_index = 0; inequality_index < Ci_rowsize; ++inequality_index) {
            a_new[x_index] += Ci[inequality_index][x_index] * extremeRay.inequality[inequality_index];
        }
    }
    double b_new = extremeRay.equality * be + extremeRay.inequality * bi;
    // store the cut
    cut_scenario.A_newRow = a_new;
    cut_scenario.b_newRow = b_new;
    return cut_scenario;
}

feasibilityCut twoStageLP_feasibilityCutGeneration_inequality(const dualMultipliers& extremeRay, const std::vector<std::vector<double>>& Ci, const std::vector<double>& bi) {
    feasibilityCut cut_scenario;
    // size of parameters
    long Ci_rowsize = Ci.size();
    long Ci_colsize = 0;
    if (Ci_rowsize > 0) {
        Ci_colsize = Ci[0].size();
    }
    long x_size = Ci_colsize; // size of x
    std::vector<double> a_new(x_size,0.0);
    for (int x_index = 0; x_index < x_size; ++x_index) {
        for (int inequality_index = 0; inequality_index < Ci_rowsize; ++inequality_index) {
            a_new[x_index] += Ci[inequality_index][x_index] * extremeRay.inequality[inequality_index];
        }
    }
    double b_new = extremeRay.inequality * bi;
    // store the cut
    cut_scenario.A_newRow = a_new;
    cut_scenario.b_newRow = b_new;
    return cut_scenario;
}

feasibilityCut twoStageLP_feasibilityCutGeneration_equality(const dualMultipliers& extremeRay, const std::vector<std::vector<double>>& Ce, const std::vector<double>& be) {
    feasibilityCut cut_scenario;
    // size of parameters
    long Ce_rowsize = Ce.size();
    long Ce_colsize = 0;
    if (Ce_rowsize > 0) {
        Ce_colsize = Ce[0].size();
    }
    long x_size = Ce_colsize; // size of x
    std::vector<double> a_new(x_size,0.0);
    for (int x_index = 0; x_index < x_size; ++x_index) {
        for (int equality_index = 0; equality_index < Ce_rowsize; ++equality_index) {
            a_new[x_index] += Ce[equality_index][x_index] * extremeRay.equality[equality_index];
        }
    }
    double b_new = extremeRay.equality * be;
    // store the cut
    cut_scenario.A_newRow = a_new;
    cut_scenario.b_newRow = b_new;
    return cut_scenario;
}

// Dx finer estimate
double DxEstimate_finer(const std::vector<double>& x, const oneStageParameters& parameters, double lowerBound, double upperBound) {
    double Dx = 0;
    // size
    long A_rowsize = parameters.A.size();
    long x_size = x.size();
    // project x onto the feasible region
    std::vector<double> x_projected = LP_projection(x, parameters, A_rowsize, x_size);
    // create an quadratic optimization problem to estimate Dx
    IloEnv env;
    IloModel mod(env);
    IloNumVarArray x_temp(env,x_size,lowerBound,upperBound,ILOFLOAT);
    IloExpr expr_obj(env);
    for (int x_index = 0; x_index < x_size; ++x_index) {
        expr_obj += x_temp[x_index] * x_temp[x_index] - 2 * x_projected[x_index] * x_temp[x_index];
    }
    IloObjective obj = IloMaximize(env,expr_obj);
    mod.add(obj);
    // constraints
    // constraints
    for (int row_index = 0; row_index < A_rowsize; ++row_index) {
        IloExpr expr_constrant(env);
        for (auto a_iterator : parameters.A[row_index]) {
            int x_index = a_iterator.first;
            double coefficient = a_iterator.second;
            //std::cout << "index of x: " << x_index << std::endl;
            expr_constrant += coefficient * x_temp[x_index];
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
    // create solver environment
    IloCplex cplex(env);
    cplex.extract(mod);
    cplex.setParam(IloCplex::OptimalityTarget, 2);
    cplex.solve();
    // calculate Dx
    for (int x_index = 0; x_index < x_size; ++x_index) {
        Dx += (cplex.getValue(x_temp[x_index])- x_projected[x_index]) * (cplex.getValue(x_temp[x_index])- x_projected[x_index]);
    }
    Dx = sqrt(Dx); // take the square root
    env.end();
    return Dx;
}

// operator overloading
double operator*(const std::vector<double>& a, const std::vector<double>& b) {
    double product = 0;
    long a_size = a.size();
    long b_size = b.size();
    if (a_size != b_size) { // if the sizes of two vectors are not equal
        std::cout << "Warning: the sizes of two vectors are not equal. Return -1" << std::endl;
        std::cout << "Warning: First vector is (";
        for (int a_index = 0; a_index < a_size; ++a_index) {
            std::cout << a[a_index] << " ";
        }
        std::cout << ")" << std::endl;
        std::cout << "Warning: Second vector is (";
        for (int b_index = 0; b_index < b_size; ++b_index) {
            std::cout << a[b_index] << " ";
        }
        std::cout << ")" << std::endl;
        return -1;
    }
    for (int index = 0; index < a_size; ++index) {
        product += a[index] * b[index];
    }
    return product;
}

// other supplemental functions
double max(const double& x, const double& y){
    if (x > y){
        return x;
    }
    else {
        return y;
    }
}

int max(const int& x, const int& y){
    if (x > y){
        return x;
    }
    else {
        return y;
    }
}

long max(const long& x, const long& y){
    if (x > y){
        return x;
    }
    else {
        return y;
    }
}
