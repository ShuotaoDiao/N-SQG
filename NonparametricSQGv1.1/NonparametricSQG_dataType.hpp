//
//  NonparametricSQG_dataType.hpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/23/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#ifndef NonparametricSQG_dataType_hpp
#define NonparametricSQG_dataType_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <ilcplex/ilocplex.h>
#include <ctime>
#include <stdlib.h>
#include <cassert>
#include <unordered_map>

// Target: nonparametricSQG
// data structures used in nonparametricSQG solver
// uppercase letter stands for matrix, lowercase letter stands for vector or scalar
// data structure for dual multipliers
struct dualMultipliers {
    std::vector<double> equality;
    std::vector<double> inequality;
    bool feasible_flag;
};
// data strcuture for the feasiblity cut
struct feasibilityCut {
    std::vector<double> A_newRow;
    double b_newRow;
};

// ****************************************************
// Target: nonparametricSQG and ioNonparametricModel
// data structure used in model setup
// data structure for the parameters in two stage linear programming
struct twoStageParameters {
    // first stage
    std::vector<double> c;
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    // second stage
    std::vector<double> d;
    std::vector<std::vector<double>> De;
    std::vector<std::vector<double>> Ce;
    std::vector<std::vector<double>> Di;
    std::vector<std::vector<double>> Ci;
};
// data structure for the parameters in one stage piecewise linear programming (hash map)
struct oneStageParameters {
    std::vector<std::string> x;
    std::vector<std::unordered_map<int,double>> A;
    std::unordered_map<int,double> b;
};

// ****************************************************
// Target: ioNonparametricDB
// data structures for the dataPoint
struct dataPoint { // definition of dataPoint
    std::vector<double> predictor;
    std::vector<double> response;
    double weight;
    // default constructor
    //dataPoint();
    // copy constructor
    //dataPoint(const dataPoint& targetPoint);
    // assignment
    //dataPoint operator=(const dataPoint& targetPoint);
};

// ****************************************************
// Target:ioStochastic
// data structure
struct randomVector {
    std::vector<double> component; // all the entries of a random vector
    std::vector<int> randomIndices; // indices of random entries
};
struct randomVectorV2 {
    std::unordered_map<int, double> component; // all the nonzero and random entries of a random vector
    std::vector<int> randomIndices; // indices of random entries
    bool flag_random = false; // flag which tells whether this vector is random
};
struct randomScalar {
    double component = 0;
    bool flag_random = false; // flag which tells whether this scalar is random
};

// vectors on the right hand side of second stage problem
struct secondStageRHS {
    randomVector be;
    randomVector bi;
};
// database of vectors on the right hand side of second stage
struct secondStageRHSDB {
    std::vector<std::vector<std::vector<double>>> be_database;
    std::vector<std::vector<std::vector<double>>> bi_database;
    std::vector<std::vector<double>> weight_database;
};

// stochastic piecewise linear objective function
// currently there are only two pieces in the max function
struct oneStagePiecewiseLinearObjective {
    randomVectorV2 c1;
    randomScalar e1;
    randomVectorV2 c2;
    randomScalar e2;
    randomVectorV2 c3;
    randomScalar e3;
    double cPiecewise = 0; // coefficient before the picewise linear function
};
// database for storing the parameters in the objective
struct oneStagePiecewiseLinearObjectiveDB {
    std::vector<std::vector<std::unordered_map<int, double>>> c1;
    std::vector<std::vector<double>> e1;
    std::vector<std::vector<std::unordered_map<int, double>>> c2;
    std::vector<std::vector<double>> e2;
    std::vector<std::vector<std::unordered_map<int, double>>> c3;
    std::vector<std::vector<double>> e3;
    std::vector<std::vector<double>> weight;
    double cPiecewise = 0; // coefficient before the piecewise linear function
};
// ****************************************************
// Target: twoStageLP_NSQG
struct validationResult {
    double mean;
    double variance;
    const double alpha = 95;
    const double Zalpha = 1.96;
    double CI_lower;
    double CI_upper;
    int num_dataPoint;
};

#endif /* NonparametricSQG_dataType_hpp */
