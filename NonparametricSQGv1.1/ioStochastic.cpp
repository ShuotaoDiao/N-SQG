//
//  ioStochastic.cpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/23/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#include "ioStochastic.hpp"

// input funtions for random vectors
secondStageRHS readStochastic(const std::string& stochasticPath) {
    secondStageRHS bRHS; // be and bi on the right hand side of constraints in the second stage problem
    const std::string nameBeginSto("<sto>");
    const std::string nameEndSto("</sto>");
    const std::string nameBeginParameter_be("<be>");
    const std::string nameEndParameter_be("</be>");
    const std::string nameBeginParameter_bi("<bi>");
    const std::string nameEndParameter_bi("</bi>");
    std::string readCondition("null"); // current condition of reading
    const char* stochasticPathConst = stochasticPath.c_str();
    std::ifstream readFile(stochasticPathConst);
    int be_curIndex = -1;
    int bi_curIndex = -1;
    if (readFile.is_open()) {
        std::string line1;
        while (getline(readFile,line1)) {
            std::stringstream ss1(line1);
            if (nameBeginSto.compare(line1) != 0 && nameEndSto.compare(line1) != 0) {
                if (nameBeginParameter_be.compare(line1) == 0) {
                    readCondition = "be"; // start reading be
                }
                else if (nameBeginParameter_bi.compare(line1) == 0) {
                    readCondition = "bi"; // start reading bi
                }
                else if (nameEndParameter_be.compare(line1) == 0) {
                    readCondition = "null"; // end reading be
                }
                else if (nameEndParameter_bi.compare(line1) == 0) {
                    readCondition = "null"; // end reading bi
                }
                else {
                    if (readCondition.compare("be") == 0) { // continue reading be
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            while (getline(ss2, line1, ',')) {
                                be_curIndex++;
                                if (line1.compare("*") == 0) {
                                    bRHS.be.randomIndices.push_back(be_curIndex); // store the index of random entry
                                    bRHS.be.component.push_back(0);
                                }
                                else { // constant entry
                                    double value;
                                    std::stringstream ss3(line1);
                                    ss3 >> value;
                                    bRHS.be.component.push_back(value);
                                }
                            }
                        }
                    }
                    else if (readCondition.compare("bi") == 0) { // continue reading bi
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            while (getline(ss2, line1, ',')) {
                                bi_curIndex++;
                                if (line1.compare("*") == 0) {
                                    bRHS.bi.randomIndices.push_back(bi_curIndex); // store the index of random entry
                                    bRHS.bi.component.push_back(0);
                                }
                                else { // constant entry
                                    double value;
                                    std::stringstream ss3(line1);
                                    ss3 >> value;
                                    bRHS.bi.component.push_back(value);
                                }
                            }
                        }
                    } // end if
                } // end if
            }// end if
        } // end while
    }
    readFile.close();
    return bRHS;
} // end readStochastic

// input function for reading random coefficients in the objective function
oneStagePiecewiseLinearObjective readStochasticOneStagePiecewiseLP(const std::string& stochasticPath) {
    oneStagePiecewiseLinearObjective parameters_objective;
    // constant parameters indicating the components of the model
    const std::string nameBeginSto("<sto>");
    const std::string nameEndSto("</sto>");
    const std::string nameBeginParameter_c1("<c1>");
    const std::string nameEndParameter_c1("</c1>");
    const std::string nameBeginParameter_e1("<e1>");
    const std::string nameEndParameter_e1("</e1>");
    const std::string nameBeginParameter_cPiecewise("<cPiecewise>");
    const std::string nameEndParameter_cPiecewise("</cPiecewise>");
    const std::string nameBeginParameter_c2("<c2>");
    const std::string nameEndParameter_c2("</c2>");
    const std::string nameBeginParameter_e2("<e2>");
    const std::string nameEndParameter_e2("</e2>");
    const std::string nameBeginParameter_c3("<c3>");
    const std::string nameEndParameter_c3("</c3>");
    const std::string nameBeginParameter_e3("<e3>");
    const std::string nameEndParameter_e3("</e3>");
    // input file setup
    std::string readCondition("null"); // current condition of reading
    const char* stochasticPathConst = stochasticPath.c_str();
    std::ifstream readFile(stochasticPathConst);
    // read file
    if (readFile.is_open()) {
        std::string line1;
        while (getline(readFile, line1)) {
            std::stringstream ss1(line1); // convert a string into a stream
            if (nameBeginSto.compare(line1) != 0 && nameEndSto.compare(line1) != 0) {
                if (nameBeginParameter_c1.compare(line1) == 0) { // start reading c1
                    readCondition = "c1";
                }
                else if (nameEndParameter_c1.compare(line1) == 0){ // end reading c1
                    readCondition = "null";
                }
                else if (nameBeginParameter_e1.compare(line1) == 0) { // start reading e1
                    readCondition = "e1";
                }
                else if (nameEndParameter_e1.compare(line1) == 0) { // end reading e1
                    readCondition = "null";
                }
                else if (nameBeginParameter_cPiecewise.compare(line1) == 0) { // start reading cPiecewise
                    readCondition = "cPiecewise";
                }
                else if (nameEndParameter_cPiecewise.compare(line1) == 0) { // end reading cPiecewise
                    readCondition = "null";
                }
                else if (nameBeginParameter_c2.compare(line1) == 0) { // start reading c2
                    readCondition = "c2";
                }
                else if (nameEndParameter_c2.compare(line1) == 0) { // end reading c2
                    readCondition = "null";
                }
                else if (nameBeginParameter_e2.compare(line1) == 0) { // start reading e2
                    readCondition = "e2";
                }
                else if (nameEndParameter_e2.compare(line1) == 0) { // end reading e2
                    readCondition = "null";
                }
                else if (nameBeginParameter_c3.compare(line1) == 0) { // start reading c3
                    readCondition = "c3";
                }
                else if (nameEndParameter_c3.compare(line1) == 0) { // end reading c3
                    readCondition = "c3";
                }
                else if (nameBeginParameter_e3.compare(line1) == 0) { // begin reading e3
                    readCondition = "e3";
                }
                else if (nameEndParameter_e3.compare(line1) == 0) { // end readning e3
                    readCondition = "null";
                }
                else {
                    if (readCondition.compare("c1") == 0) { // continue reading c1
                        while (getline(ss1, line1, ';')) {
                            //std::cout << "Debug: " << line1 << std::endl;
                            std::stringstream ss2(line1); // convert a string into a stream
                            int c1_index = 0;
                            double c1_value = 0;
                            int index = 0;
                            while (getline(ss2, line1, ':')) {
                                //std::cout << "Debug: " << line1 << std::endl;
                                std::stringstream ss3(line1);
                                if (index == 0) { // read c1 index
                                    ss3 >> c1_index;
                                }
                                else { // read value at current c1 index
                                    if (line1.compare("*") == 0) { // if it is random
                                        c1_value = 0;
                                        parameters_objective.c1.randomIndices.push_back(c1_index); // store the index of random entry
                                        if (parameters_objective.c1.flag_random == false) {
                                            parameters_objective.c1.flag_random = true; // tell the vetor is random
                                        }
                                    }
                                    else { // if it is constant
                                        ss3 >> c1_value;
                                    }
                                }
                                index++;
                            }
                            // store
                            //std::cout << c1_index << " " << c1_value << std::endl;
                            parameters_objective.c1.component[c1_index] = c1_value;
                        } // end while
                    } // end if
                    else if (readCondition.compare("e1") == 0) { // continue reading e1
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1); // convert a string into a stream
                            double e1_value = 0;
                            if (line1.compare("*") == 0) { // e1 is random
                                parameters_objective.e1.component = e1_value;
                                parameters_objective.e1.flag_random = true;
                            }
                            else { // e1 is a constant
                                ss2 >> e1_value;
                                parameters_objective.e1.component = e1_value;
                            }
                        }
                    }
                    else if (readCondition.compare("cPiecewise") == 0) { // continue reading cPiecewise
                        std::stringstream ss2(line1); // convert a string into a stream
                        double cPicewise_value = 0;
                        ss2 >> cPicewise_value;
                        parameters_objective.cPiecewise = cPicewise_value;
                    }
                    else if (readCondition.compare("c2") == 0) { // continue reading c2
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1); // convert a string into a stream
                            int c2_index = 0;
                            double c2_value = 0;
                            int index = 0;
                            while (getline(ss2, line1, ':')) {
                                std::stringstream ss3(line1); // convert a string into a stream
                                if (index == 0) { // read c2 index
                                    ss3 >> c2_index;
                                }
                                else { // read value at current c2 index
                                    if (line1.compare("*") == 0) { // if it is random
                                        c2_value = 0;
                                        parameters_objective.c2.randomIndices.push_back(c2_index); // store the index of random entry
                                        if (parameters_objective.c2.flag_random == false) {
                                            parameters_objective.c2.flag_random = true; // tell the vetor is random
                                        }
                                    }
                                    else { // it is constant
                                        ss3 >> c2_value;
                                    }
                                }
                                index++;
                            }
                            // store
                            parameters_objective.c2.component[c2_index] = c2_value;
                        }
                    }
                    else if (readCondition.compare("e2") == 0) { // continue reading e2
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1); // convert a string into a stream
                            double e2_value = 0;
                            if (line1.compare("*") == 0) { // e2 is random
                                parameters_objective.e2.component = e2_value;
                                parameters_objective.e2.flag_random = true;
                            }
                            else { // e2 is a constant
                                ss2 >> e2_value;
                                parameters_objective.e2.component = e2_value;
                            }
                        }
                    }
                    else if (readCondition.compare("c3") == 0) { // continue reading c3
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1); // convert a string into a stream
                            int c3_index = 0;
                            double c3_value = 0;
                            int index = 0;
                            while (getline(ss2, line1, ':')) {
                                std::stringstream ss3(line1); // convert a string into a stream
                                if (index == 0) { // read c3 index
                                    ss3 >> c3_index;
                                }
                                else { // read value at current c3 index
                                    if (line1.compare("*") == 0) { // if it is random
                                        c3_value = 0;
                                        parameters_objective.c3.randomIndices.push_back(c3_index); // store the index of random entry
                                        if (parameters_objective.c3.flag_random == false) {
                                            parameters_objective.c3.flag_random = true; // tell the vetor is random
                                        }
                                    }
                                    else {
                                        ss3 >> c3_value;
                                    }
                                }
                                index++;
                            }
                            // store
                            parameters_objective.c3.component[c3_index] = c3_value;
                        }
                    }
                    else if (readCondition.compare("e3") == 0) { // continue reading e3
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1); // convert a string into a stream
                            double e3_value = 0;
                            if (line1.compare("*") == 0) { // e3 is random
                                parameters_objective.e3.component = e3_value;
                                parameters_objective.e3.flag_random = true;
                            }
                            else { // e3 is a constant
                                ss2 >> e3_value;
                                parameters_objective.e3.component = e3_value;
                            }
                        }
                    }
                } // end if
            } // end if
        } // end while
    } // end if
    else {
        std::cout << "ERROR: sto.txt is not found!" << std::endl;
    }
    return parameters_objective;
}

// functions for merging nonparametericDB and randomVector
secondStageRHSDB mergeDB_randomVector(const std::vector<std::vector<dataPoint>>& be_DB, const std::vector<std::vector<dataPoint>>& bi_DB, const secondStageRHS& bRHS) {
    secondStageRHSDB bRHSDB;
    // number of datasets in a database
    long be_numOfDataset = be_DB.size();
    long bi_numOfDataset = bi_DB.size();
    std::cout << "be_numOfDataset: " << be_numOfDataset << std::endl;
    std::cout << "bi_numOfDataset: " << bi_numOfDataset << std::endl;
    // generate be
    for (int be_dataset_index = 0; be_dataset_index < be_numOfDataset; ++be_dataset_index) {
        std::vector<std::vector<double>> be_dataset_temp;
        std::vector<double> weight_dataset_temp;
        // number of datapoints in the current dataset
        long be_numOfDatapoint = be_DB[be_dataset_index].size();
        for (int be_datapoint_index = 0; be_datapoint_index < be_numOfDatapoint; ++be_datapoint_index) {
            // obtain weight of current datapoint
            weight_dataset_temp.push_back(be_DB[be_dataset_index][be_datapoint_index].weight);
            // number of entries in a datapoint
            long be_numOfEntry = bRHS.be.component.size();
            // initializaion
            std::vector<double> be_temp(be_numOfEntry,0.0);
            int be_random_index = 0;
            long be_numOfRandom = bRHS.be.randomIndices.size();
            for (int entry_index = 0; entry_index < be_numOfEntry; ++entry_index) {
                if (be_random_index < be_numOfRandom) {
                    if (entry_index == bRHS.be.randomIndices[be_random_index]) { // if current entry is random
                        be_temp[entry_index] = be_DB[be_dataset_index][be_datapoint_index].response[be_random_index];
                        be_random_index++;
                    }
                    else {
                        be_temp[entry_index] = bRHS.be.component[entry_index];
                    }
                }
                else {
                    be_temp[entry_index] = bRHS.be.component[entry_index];
                } // end if
            } // end for loop
            be_dataset_temp.push_back(be_temp);
        } // end for loop
        bRHSDB.be_database.push_back(be_dataset_temp);
        bRHSDB.weight_database.push_back(weight_dataset_temp);
    }
    if (be_numOfDataset < 1) {
        // generate bi
        for (int bi_dataset_index = 0; bi_dataset_index < bi_numOfDataset; ++bi_dataset_index) {
            std::cout << "bi_dataste_index: " << bi_dataset_index << std::endl;
            std::vector<std::vector<double>> bi_dataset_temp;
            std::vector<double> weight_dataset_temp;
            // number of datapoints in the current dataset
            long bi_numOfDatapoint = bi_DB[bi_dataset_index].size();
            for (int bi_datapoint_index = 0; bi_datapoint_index < bi_numOfDatapoint; ++bi_datapoint_index) {
                std::cout << "bi_datapoint_index: " << bi_datapoint_index << std::endl;
                // obtain weight of current datapoint
                weight_dataset_temp.push_back(bi_DB[bi_dataset_index][bi_datapoint_index].weight);
                // number of entries in a datapoint
                long bi_numOfEntry = bRHS.bi.component.size();
                // initializaion
                std::vector<double> bi_temp(bi_numOfEntry,0.0);
                int bi_random_index = 0;
                long bi_numOfRandom = bRHS.bi.randomIndices.size();
                for (int entry_index = 0; entry_index < bi_numOfEntry; ++entry_index) {
                    //std::cout << "entry_index: " << entry_index << std::endl;
                    if (bi_random_index < bi_numOfRandom) {
                        if (entry_index == bRHS.bi.randomIndices[bi_random_index]) { // if current entry is random
                            bi_temp[entry_index] = bi_DB[bi_dataset_index][bi_datapoint_index].response[bi_random_index];
                            bi_random_index++;
                        }
                        else {
                            bi_temp[entry_index] = bRHS.bi.component[entry_index];
                        }
                    }
                    else {
                        bi_temp[entry_index] = bRHS.bi.component[entry_index];
                    } // end if
                } // end for loop
                bi_dataset_temp.push_back(bi_temp);
            } // end for loop
            bRHSDB.bi_database.push_back(bi_dataset_temp);
            bRHSDB.weight_database.push_back(weight_dataset_temp);
        }
    }
    else {
        // generate bi
        for (int bi_dataset_index = 0; bi_dataset_index < bi_numOfDataset; ++bi_dataset_index) {
            std::vector<std::vector<double>> bi_dataset_temp;
            // number of datapoints in the current dataset
            long bi_numOfDatapoint = bi_DB[bi_dataset_index].size();
            for (int bi_datapoint_index = 0; bi_datapoint_index < bi_numOfDatapoint; ++bi_datapoint_index) {
                // number of entries in a datapoint
                long bi_numOfEntry = bRHS.bi.component.size();
                // initializaion
                std::vector<double> bi_temp(bi_numOfEntry,0.0);
                int bi_random_index = 0;
                long bi_numOfRandom = bRHS.be.randomIndices.size();
                for (int entry_index = 0; entry_index < bi_numOfEntry; ++entry_index) {
                    if (bi_random_index < bi_numOfRandom) {
                        if (entry_index == bRHS.bi.randomIndices[bi_random_index]) { // if current entry is random
                            bi_temp[entry_index] = bi_DB[bi_dataset_index][bi_datapoint_index].response[bi_random_index];
                            bi_random_index++;
                        }
                        else {
                            bi_temp[entry_index] = bRHS.bi.component[entry_index];
                        }
                    }
                    else {
                        bi_temp[entry_index] = bRHS.bi.component[entry_index];
                    } // end if
                } // end for loop
                bi_dataset_temp.push_back(bi_temp);
            } // end for loop
            bRHSDB.bi_database.push_back(bi_dataset_temp);
        }
    }
    return bRHSDB;
} // end mergeDB_randomVector

// functions for mergeing nonparametricDB and random objetive function
oneStagePiecewiseLinearObjectiveDB mergeDB_randomObjective(const std::vector<std::vector<dataPoint>>& DB, oneStagePiecewiseLinearObjective parameters_objective, long x_size) {
    oneStagePiecewiseLinearObjectiveDB objectiveDB;
    // number of dataset
    long num_datasets = DB.size();
    // number of random entries in each vector
    long num_random_c1 = parameters_objective.c1.randomIndices.size();
    long num_random_c2 = parameters_objective.c2.randomIndices.size();
    long num_random_c3 = parameters_objective.c3.randomIndices.size();
    // store the coefficient before the piecewise linear function
    objectiveDB.cPiecewise = parameters_objective.cPiecewise;
    // transform dataset based on the locations of random entries
    for (int dataset_index = 0; dataset_index < num_datasets; ++dataset_index) {
        //std::cout << "Debugging dataset_index : " << dataset_index << std::endl;
        // dataset for c1
        std::vector<std::unordered_map<int, double>> c1_dataset;
        // dataset for e1
        std::vector<double> e1_dataset;
        // dataset for c2
        std::vector<std::unordered_map<int, double>> c2_dataset;
        // dataset for e2
        std::vector<double> e2_dataset;
        // dataset for c3
        std::vector<std::unordered_map<int, double>> c3_dataset;
        // dataset for e3
        std::vector<double> e3_dataset;
        // dataset for weights
        std::vector<double> weight_dataset;
        // number of scenarios
        long num_scenarios = DB[dataset_index].size();
        // loop through all the scenarios (data points)
        for (int scenario_index = 0; scenario_index < num_scenarios; ++scenario_index) {
            weight_dataset.push_back(DB[dataset_index][scenario_index].weight); // add weight
            // c1 c2 c3
            std::unordered_map<int, double> c1_datapoint;
            std::unordered_map<int, double> c2_datapoint;
            std::unordered_map<int, double> c3_datapoint;
            int base_index = 0;
            int random_index = 0;
            // store c1
            for (int x_index = 0; x_index < x_size; ++x_index) {
                if (parameters_objective.c1.component.find(x_index) != parameters_objective.c1.component.end()) { // coefficient for x[x_index] is nonzero or random
                    if (random_index < num_random_c1) { // the coefficient is random
                        if (x_index == parameters_objective.c1.randomIndices[random_index]) {
                            c1_datapoint[x_index] = DB[dataset_index][scenario_index].response[base_index + random_index];
                            random_index++;
                        }
                        else { // the coefficient is a constant
                            c1_datapoint[x_index] = parameters_objective.c1.component[x_index];
                        }
                    }
                    else { // the coefficient is a constant
                        c1_datapoint[x_index] = parameters_objective.c1.component[x_index];
                    }
                }
            } // end for
            // refresh random index and update base_index
            base_index = random_index;
            random_index = 0;
            if (parameters_objective.e1.flag_random == true) { // if e1 is random
                e1_dataset.push_back(DB[dataset_index][scenario_index].response[base_index + 1]);
                base_index++;
            }
            else { // if e1 is constant
                e1_dataset.push_back(parameters_objective.e1.component);
            } // end if
            // store c2
            for (int x_index = 0; x_index < x_size; ++x_index) {
                if (parameters_objective.c2.component.find(x_index) != parameters_objective.c2.component.end()) { // coefficient for x[x_index] is nonzero or random
                    if (random_index < num_random_c2) { // the coefficient is random
                        if (x_index == parameters_objective.c2.randomIndices[random_index]) {
                            c2_datapoint[x_index] = DB[dataset_index][scenario_index].response[base_index + random_index]; // store the coefficient
                            random_index++;
                        }
                        else { // the coefficient is a constant
                            c2_datapoint[x_index] = parameters_objective.c2.component[x_index]; // store the coefficient
                        }
                    }
                    else { // the coefficient is a constant
                        c2_datapoint[x_index] = parameters_objective.c2.component[x_index]; // store the coefficient
                    }
                }
            } // end for
            // refresh random index and update base_index
            base_index = random_index;
            random_index = 0;
            if (parameters_objective.e2.flag_random == true) { // if e2 is random
                e2_dataset.push_back(DB[dataset_index][scenario_index].response[base_index + 1]);
                base_index++;
            }
            else { // if e1 is constant
                e2_dataset.push_back(parameters_objective.e2.component);
            } // end if
            // store c3
            for (int x_index = 0; x_index < x_size; ++x_index) {
                if (parameters_objective.c3.component.find(x_index) != parameters_objective.c3.component.end()) { // coefficient for x[x_index] is nonzero or random
                    if (random_index < num_random_c3) { // the coefficient is random
                        if (x_index == parameters_objective.c3.randomIndices[random_index]) {
                            c3_datapoint[x_index] = DB[dataset_index][scenario_index].response[base_index + random_index]; // store the coefficient
                            random_index++;
                        }
                        else { // the coefficient is a constant
                            c3_datapoint[x_index] = parameters_objective.c3.component[x_index]; // store the coefficient
                        }
                    }
                    else { // the coefficient is a constant
                        c3_datapoint[x_index] = parameters_objective.c3.component[x_index]; // store the coefficient
                    }
                }
            } // end for
            // refresh random index and update base_index
            base_index = random_index;
            random_index = 0;
            if (parameters_objective.e3.flag_random == true) { // if e3 is random
                e3_dataset.push_back(DB[dataset_index][scenario_index].response[base_index + 1]);
                base_index++;
            }
            else { // if e1 is constant
                e3_dataset.push_back(parameters_objective.e3.component);
            } // end if
            // store the datapoints into dataset
            c1_dataset.push_back(c1_datapoint);
            c2_dataset.push_back(c2_datapoint);
            c3_dataset.push_back(c3_datapoint);
        }
        // store the datasets into objectiveDB
        objectiveDB.c1.push_back(c1_dataset);
        objectiveDB.e1.push_back(e1_dataset);
        objectiveDB.c2.push_back(c2_dataset);
        objectiveDB.e2.push_back(e2_dataset);
        objectiveDB.c3.push_back(c3_dataset);
        objectiveDB.e3.push_back(e3_dataset);
        objectiveDB.weight.push_back(weight_dataset);
    }
    return objectiveDB;
} // end mergeDB_randomObjective


// output functions
void printStochastic(const secondStageRHS& bRHS) {
    // sizes of be and bi
    long be_size = bRHS.be.component.size();
    long bi_size = bRHS.bi.component.size();
    // number of random entries
    long beRandom_size = bRHS.be.randomIndices.size();
    long biRandom_size = bRHS.bi.randomIndices.size();
    std::cout << "biRandom_size : " << biRandom_size << std::endl;
    // pirnt be
    if (be_size > 0) {
        std::cout << "be: ";
        int random_curIndex = 0;
        for (int index = 0; index < be_size - 1; ++index) {
            if (random_curIndex < beRandom_size) {
                if (index == bRHS.be.randomIndices[random_curIndex]) { // current entry is random
                    random_curIndex++;
                    std::cout << "*,";
                }
                else {
                    std::cout << bRHS.be.component[index] << ",";
                }
            }
            else {
                std::cout << bRHS.be.component[index] << ",";
            }
        } // end for loop
        int index = be_size - 1;
        if (random_curIndex < beRandom_size) {
            if (index == bRHS.be.randomIndices[random_curIndex]) { // current entry is random
                random_curIndex++;
                std::cout << "*,";
            }
            else {
                std::cout << bRHS.be.component[index] << std::endl;
            }
        }
        else {
            std::cout << bRHS.be.component[index] << std::endl;
        }
    }
    else {
        std::cout << "be: NULL" << std::endl;
    } // end if
    // print bi
    if (bi_size > 0) {
        std::cout << "bi: ";
        int random_curIndex = 0;
        for (int index = 0; index < bi_size - 1; ++index) {
            if (random_curIndex < biRandom_size) {
                if (index == bRHS.bi.randomIndices[random_curIndex]) { // current entry is random
                    random_curIndex++;
                    std::cout << "*,";
                }
                else {
                    std::cout << bRHS.bi.component[index] << ",";
                }
            }
            else {
                std::cout << bRHS.bi.component[index] << ",";
            }
        } // end for loop
        int index = bi_size - 1;
        if (random_curIndex < biRandom_size) {
            if (index == bRHS.bi.randomIndices[random_curIndex]) { // current entry is random
                random_curIndex++;
                std::cout << "*,";
            }
            else {
                std::cout << bRHS.bi.component[index] << std::endl;
            }
        }
        else {
            std::cout << bRHS.bi.component[index] << std::endl;
        }
    }
    else {
        std::cout << "bi: NULL" << std::endl;
    }
}// end printStochastic

// print the paramerts of the piecewise linear objective function
void printStochasticOneStagePiecewiseLP(oneStagePiecewiseLinearObjective parameters_objective, oneStageParameters parameters) {
    // size
    long x_size = parameters.x.size(); // size of x
    // print c1
    std::cout << ">>>c1<<<" << std::endl;
    for (int x_index = 0; x_index < x_size; ++x_index) {
        int random_index = 0;
        long randomIndices_size = parameters_objective.c1.randomIndices.size();
        if (parameters_objective.c1.component.find(x_index) != parameters_objective.c1.component.end()) { // coefficient for x[x_index] is nonzero or random
            if (random_index < randomIndices_size) { // the coefficient is random
                if (x_index == parameters_objective.c1.randomIndices[random_index]) {
                    std::cout << "Position: " << x_index << " Coefficient: random" << std::endl;
                    random_index++;
                }
                else { // the coefficient is a constant
                    std::cout << "Position: " << x_index << " Coefficient: " << parameters_objective.c1.component[x_index] << std::endl;
                }
            }
            else { // the coefficient is a constant
                std::cout << "Position: " << x_index << " Coefficient: " << parameters_objective.c1.component[x_index] << std::endl;
            }
        }
        else {
            std::cout << "Position: " << x_index << " Coefficient: 0" << std::endl;
        }
    }
    // print e1
    std::cout << ">>>e1<<<" << std::endl;
    if (parameters_objective.e1.flag_random == true) { // if it is random
        std::cout << "Scalar: random" << std::endl;
    }
    else {
        std::cout << "Scalar: " << parameters_objective.e1.component << std::endl;
    }
    // print c2
    std::cout << ">>>c2<<<" << std::endl;
    for (int x_index = 0; x_index < x_size; ++x_index) {
        int random_index = 0;
        long randomIndices_size = parameters_objective.c2.randomIndices.size();
        if (parameters_objective.c2.component.find(x_index) != parameters_objective.c2.component.end()) {// the coefficient for x[x_index] is nonzero
            if (random_index < randomIndices_size) { // the coefficient is random
                if (x_index == parameters_objective.c2.randomIndices[random_index]) {
                    std::cout << "Position: " << x_index << " Coefficient: random" << std::endl;
                    random_index++;
                }
                else { // the coefficient is a constant
                    std::cout << "Position: " << x_index << " Coefficient: " << parameters_objective.c2.component[x_index] << std::endl;
                }
            }
            else { // the coefficient is a constant
                std::cout << "Position: " << x_index << " Coefficient: " << parameters_objective.c2.component[x_index] << std::endl;
            }
        }
        else {
            std::cout << "Position: " << x_index << " Coefficient: 0" << std::endl;
        }
    }
    // print e2
    std::cout << ">>>e2<<<" << std::endl;
    if (parameters_objective.e2.flag_random == true) { // if it is random
        std::cout << "Scalar: random" << std::endl;
    }
    else {
        std::cout << "Scalar: " << parameters_objective.e2.component << std::endl;
    }
    // print c3
    std::cout << ">>>c3<<<" << std::endl;
    for (int x_index = 0; x_index < x_size; ++x_index) {
        int random_index = 0;
        long randomIndices_size = parameters_objective.c3.randomIndices.size();
        if (parameters_objective.c3.component.find(x_index) != parameters_objective.c3.component.end()) {// the coefficient for x[x_index] is nonzero
            if (random_index < randomIndices_size) { // the coefficient is random
                if (x_index == parameters_objective.c3.randomIndices[random_index]) {
                    std::cout << "Position: " << x_index << " Coefficient: random" << std::endl;
                    random_index++;
                }
                else { // the coefficient is a constant
                    std::cout << "Position: " << x_index << " Coefficient: " << parameters_objective.c3.component[x_index] << std::endl;
                }
            }
            else { // the coefficient is a constant
                std::cout << "Position: " << x_index << " Coefficient: " << parameters_objective.c3.component[x_index] << std::endl;
            }
        }
        else {
            std::cout << "Position: " << x_index << " Coefficient: 0" << std::endl;
        }
    }
    // print e3
    std::cout << ">>>e3<<<" << std::endl;
    if (parameters_objective.e3.flag_random == true) { // if it is random
        std::cout << "Scalar: random" << std::endl;
    }
    else {
        std::cout << "Scalar: " << parameters_objective.e3.component << std::endl;
    }
}

