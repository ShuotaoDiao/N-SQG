//
//  ioNonparametricModel.cpp
//  NonparametricSQGv1.1
//
//  Created by Shuotao Diao on 6/23/19.
//  Copyright © 2019 Shuotao Diao. All rights reserved.
//

#include "ioNonparametricModel.hpp"

//input functions for two stage linear programming
twoStageParameters readTwoStageParameters(const std::string& parameterPath) {
    twoStageParameters parameters;
    const std::string nameBeginModel("<model:twoStageLP>");
    const std::string nameEndModel("</model:twoStageLP>");
    const std::string nameBeginParameter_c("<c>");
    const std::string nameEndParameter_c("</c>");
    const std::string nameBeginParameter_A("<A>");
    const std::string nameEndParameter_A("</A>");
    const std::string nameBeginParameter_b("<b>");
    const std::string nameEndParameter_b("</b>");
    const std::string nameBeginParameter_d("<d>");
    const std::string nameEndParameter_d("</d>");
    const std::string nameBeginParameter_Di("<Di>");
    const std::string nameEndParameter_Di("</Di>");
    const std::string nameBeginParameter_Ci("<Ci>");
    const std::string nameEndParameter_Ci("</Ci>");
    const std::string nameBeginParameter_De("<De>");
    const std::string nameEndParameter_De("</De>");
    const std::string nameBeginParameter_Ce("<Ce>");
    const std::string nameEndParameter_Ce("</Ce>");
    std::string readCondition("null");
    const char* parameterPathConst = parameterPath.c_str(); // convert a string path into a constant path
    std::ifstream readFile(parameterPathConst);
    if (readFile.is_open()) {
        std::string line1;
        while (getline(readFile,line1)) {// get the whole line
            std::cout << line1 << std::endl;
            std::stringstream ss1(line1); // convert a string into stream
            if (nameBeginModel.compare(line1) != 0 && nameEndModel.compare(line1) != 0) { // main content
                if (nameBeginParameter_c.compare(line1) == 0) { // beign reading parameter c
                    readCondition = "c";
                }
                else if (nameBeginParameter_A.compare(line1) == 0) { // begin reading parameter A
                    readCondition = "A";
                }
                else if (nameBeginParameter_b.compare(line1) == 0) { // begin reading parameter b
                    readCondition = "b";
                }
                else if (nameBeginParameter_d.compare(line1) == 0) { // begin reading parameter d
                    readCondition = "d";
                }
                else if (nameBeginParameter_Di.compare(line1) == 0) { // begin reading parameter Di
                    readCondition = "Di";
                }
                else if (nameBeginParameter_Ci.compare(line1) == 0) { // begin reading parameter Ci
                    readCondition = "Ci";
                }
                else if (nameBeginParameter_De.compare(line1) == 0) { // begin reading parameter De
                    readCondition = "De";
                }
                else if (nameBeginParameter_Ce.compare(line1) == 0) { // begin reading parameter Ce
                    readCondition = "Ce";
                }
                else if (nameEndParameter_c.compare(line1) == 0) { // end reading parameter c
                    readCondition = "null";
                }
                else if (nameEndParameter_A.compare(line1) == 0) { // end reading parameter A
                    readCondition = "null";
                }
                else if (nameEndParameter_b.compare(line1) == 0) { // end reading parameter b
                    readCondition = "null";
                }
                else if (nameEndParameter_d.compare(line1) == 0) { // end reading parameter d
                    readCondition = "null";
                }
                else if (nameEndParameter_Di.compare(line1) == 0) { // end reading parameter Di
                    readCondition = "null";
                }
                else if (nameEndParameter_Ci.compare(line1) == 0) { // end reading parameter Ci
                    readCondition = "null";
                }
                else if (nameEndParameter_De.compare(line1) == 0) { // end reading parameter De
                    readCondition = "null";
                }
                else if (nameEndParameter_Ce.compare(line1) == 0) { // end reading parameter Ce
                    readCondition = "null";
                }
                else {
                    if (readCondition.compare("c") == 0) {
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            while (getline(ss2, line1, ',')) {
                                double value;
                                std::stringstream ss3(line1);
                                ss3 >> value;
                                parameters.c.push_back(value);
                            }
                        }
                    }
                    else if (readCondition.compare("A") == 0) {
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            std::vector<double> A_row; // create one new row of A matrix
                            while (getline(ss2, line1, ',')) {
                                double value;
                                std::stringstream ss3(line1);
                                ss3 >> value;
                                A_row.push_back(value);
                            }
                            parameters.A.push_back(A_row);
                        }
                    }
                    else if (readCondition.compare("b") == 0) {
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            while (getline(ss2, line1, ',')) {
                                double value;
                                std::stringstream ss3(line1);
                                ss3 >> value;
                                parameters.b.push_back(value);
                            }
                        }
                    }
                    else if (readCondition.compare("d") == 0) {
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            while (getline(ss2, line1, ',')) {
                                double value;
                                std::stringstream ss3(line1);
                                ss3 >> value;
                                parameters.d.push_back(value);
                            }
                        }
                    }
                    else if (readCondition.compare("Di") == 0) {
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            std::vector<double> Di_row; // create one new row of Di matrix
                            while (getline(ss2, line1, ',')) {
                                double value;
                                std::stringstream ss3(line1);
                                ss3 >> value;
                                Di_row.push_back(value);
                            }
                            parameters.Di.push_back(Di_row);
                        }
                    }
                    else if (readCondition.compare("Ci") == 0) {
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            std::vector<double> Ci_row; // create one new row of Ci matrix
                            while (getline(ss2, line1, ',')) {
                                double value;
                                std::stringstream ss3(line1);
                                ss3 >> value;
                                Ci_row.push_back(value);
                            }
                            parameters.Ci.push_back(Ci_row);
                        }
                    }
                    else if (readCondition.compare("De") == 0) {
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            std::vector<double> De_row; // create one new row of De matrix
                            while (getline(ss2, line1, ',')) {
                                double value;
                                std::stringstream ss3(line1);
                                ss3 >> value;
                                De_row.push_back(value);
                            }
                            parameters.De.push_back(De_row);
                        }
                    }
                    else if (readCondition.compare("Ce") == 0) {
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1);
                            std::vector<double> Ce_row; // create one new row of Ce matrix
                            while (getline(ss2, line1, ',')) {
                                double value;
                                std::stringstream ss3(line1);
                                ss3 >> value;
                                Ce_row.push_back(value);
                            }
                            parameters.Ce.push_back(Ce_row);
                        }
                    }
                }
            }
            std::cout << "Condition : " << readCondition << std::endl;
        }
    }
    readFile.close();
    
    return parameters;
}

// input functions for one stage piecewise linear programming
oneStageParameters readOneStageParameters(const std::string& parameterPath) {
    oneStageParameters parameters;
    // constant parameters indicating the components of a model
    const std::string nameBeginModel("<model:oneStagePiecewiseLP>");
    const std::string nameEndModel("</model:oneStagePiecewiseLP>");
    const std::string nameBeginParameter_x("<DecisionVariable>");
    const std::string nameEndParameter_x("</DecisionVariable>");
    const std::string nameBeginParameter_A("<A>");
    const std::string nameEndParameter_A("</A>");
    const std::string nameBeginParameter_b("<b>");
    const std::string nameEndParameter_b("</b>");
    // input setup
    std::string readCondition("null");
    const char* parameterPathConst = parameterPath.c_str(); // convert a string path into a constant path
    std::ifstream readFile(parameterPathConst);
    if (readFile.is_open()) {
        std::string line1;
        while (getline(readFile,line1)) {// get the whole line
            std::stringstream ss1(line1); // convert a string into stream
            if (nameBeginModel.compare(line1) != 0 && nameEndModel.compare(line1) != 0) { // main content
                if (nameBeginParameter_A.compare(line1) == 0) { // start reading A
                    readCondition = "A";
                }
                else if (nameEndParameter_A.compare(line1) == 0) { // end reading A
                    readCondition = "null";
                }
                else if (nameBeginParameter_b.compare(line1) == 0) { // start reading b
                    readCondition = "b";
                }
                else if (nameEndParameter_b.compare(line1) == 0) { // end reading b
                    readCondition = "null";
                }
                else if (nameBeginParameter_x.compare(line1) == 0) { // start reading decision variable
                    readCondition = "DecisionVariable";
                }
                else if (nameEndParameter_x.compare(line1) == 0) { // end reading decision variable
                    readCondition = "null";
                }
                else {
                    if (readCondition.compare("DecisionVariable") == 0) { // continue reading decision variable
                        while (getline(ss1, line1, ';')) { // get one decision variable
                            std::stringstream ss2(line1); // convert a string into a stream
                            int index1 = 0;
                            while (getline(ss2, line1, ':')) {
                                if (index1 == 1) {
                                    parameters.x.push_back(line1);
                                }
                                index1++;
                            }
                        }
                    }
                    else if (readCondition.compare("A") == 0) { // continue reading A
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1); // convert a string into a stream
                            int index1 = 0;
                            int row_index = 0; // position
                            int col_index = 0;
                            double value = 0;
                            // read the key value pair
                            while (getline(ss2, line1, ':')) {
                                if (index1 == 0) { // the position for key
                                    std::stringstream ss3(line1); // convert a string into a stream
                                    int index2 = 0;
                                    while (getline(ss3, line1, ',')) {
                                        long line1_size = line1.size();
                                        // convert a string into an integer
                                        if (index2 == 0) { // e.g. (1
                                            for (int index3 = 1; index3 < line1_size; ++index3) {
                                                int digit = line1.at(index3) - '0';
                                                row_index = row_index * 10 + digit;
                                            }
                                        }
                                        else { // e.g. 1)
                                            for (int index3 = 0; index3 < line1_size - 1; ++index3) {
                                                int digit = line1.at(index3) - '0';
                                                col_index = col_index * 10 + digit;
                                            }
                                        } // end if
                                        index2++;
                                    }
                                }
                                else {
                                    std::stringstream ss3(line1);
                                    //std::cout << "Debug Point" << std::endl;
                                    //std::cout << line1 << std::endl;
                                    ss3 >> value;
                                }
                                index1++;
                            }
                            // store the key value pair into A
                            // number of rows in A
                            //std::cout << "row : " << row_index << std::endl;
                            //std::cout << "col : " << col_index << std::endl;
                            //std::cout << "value : " << value << std::endl;
                            long A_rowsize = parameters.A.size();
                            if (row_index < A_rowsize) {
                                parameters.A[row_index][col_index] = value;
                            }
                            else {
                                while (row_index >= A_rowsize) {
                                    std::unordered_map<int, double> A_row;
                                    parameters.A.push_back(A_row);
                                    A_rowsize++;
                                }
                                parameters.A[row_index][col_index] = value;
                            }
                        } // end while
                    }
                    else if (readCondition.compare("b") == 0) { // continue reading b
                        while (getline(ss1, line1, ';')) {
                            std::stringstream ss2(line1); // convert a string into stream
                            int b_index = 0;
                            double value = 0;
                            int index1 = 0;
                            while (getline(ss2, line1, ':')) {
                                if (index1 == 0) { // information of the position of b
                                    std::stringstream ss3(line1);
                                    ss3 >> b_index;
                                }
                                else { // value at the current position of b
                                    std::stringstream ss3(line1);
                                    ss3 >> value;
                                }
                                index1++;
                            }
                            parameters.b[b_index] = value;
                        }
                    }
                } // end if
            } // end if
        } // end while
    } // end if
    else {
        std::cout << "ERROR: model.txt is not found!" << std::endl;
    }
    return parameters;
}

void printTwoStageParameters(const twoStageParameters& parameters) {
    // obtain the sizes parameters
    long c_size = parameters.c.size();
    long A_rowsize = parameters.A.size();
    long A_colsize = 0;
    if (A_rowsize > 0) {
        A_colsize = parameters.A[0].size();
    }
    long b_size = parameters.b.size();
    long d_size = parameters.d.size();
    long Di_rowsize = parameters.Di.size();
    long Di_colsize = 0;
    if (Di_rowsize > 0) {
        Di_colsize = parameters.Di[0].size();
    }
    long Ci_rowsize = parameters.Ci.size();
    long Ci_colsize = 0;
    if (Ci_rowsize > 0) {
        Ci_colsize = parameters.Ci[0].size();
    }
    std::cout << "Ci_rowsize: " << Ci_rowsize << std::endl;
    // print c
    std::cout << "c: ";
    for (int index = 0; index < c_size; ++index) {
        std::cout << parameters.c[index] << " ";
    }
    std::cout << std::endl;
    // print A
    if (A_rowsize > 0) {
        std::cout << "A: ";
        for (int row_index = 0; row_index < A_rowsize; ++row_index) {
            for (int col_index = 0; col_index < A_colsize; ++col_index) {
                std::cout << parameters.A[row_index][col_index] << " ";
            }
            std::cout << std::endl;
        }
    }
    else {
        std::cout << "A: NULL" << std::endl;
    }
    // print b
    if (b_size > 0) {
        std::cout << "b: ";
        for (int index = 0; index < b_size; ++index) {
            std::cout << parameters.b[index] << " ";
        }
        std::cout << std::endl;
    }
    else {
        std::cout << "b: NULL" << std::endl;
    }
    // print d
    if (d_size > 0) {
        std::cout << "d: ";
        for (int index = 0; index < d_size; ++index) {
            std::cout << parameters.d[index] << " ";
        }
        std::cout << std::endl;
    }
    else {
        std::cout << "d: NULL" << std::endl;
    }
    // print Di
    if (Di_rowsize > 0) {
        std::cout << "Di: ";
        for (int row_index = 0; row_index < Di_rowsize; ++row_index) {
            for (int col_index = 0; col_index < Di_colsize; ++col_index) {
                std::cout << parameters.Di[row_index][col_index] << " ";
            }
            std::cout << std::endl;
        }
    }
    else {
        std::cout << "Di: NULL" << std::endl;
    }
    // print Ci
    if (Di_rowsize > 0) {
        std::cout << "Ci: ";
        for (int row_index = 0; row_index < Ci_rowsize; ++row_index) {
            for (int col_index = 0; col_index < Ci_colsize; ++col_index) {
                std::cout << parameters.Ci[row_index][col_index] << " ";
            }
            std::cout << std::endl;
        }
    }
    else {
        std::cout << "Ci: NULL" << std::endl;
    }
}

// print one stage parameters (A and b)
void printOneStageParameters(oneStageParameters parameters) {
    // size of A and b
    long x_size = parameters.x.size(); // number of decision variables
    long A_rowsize = parameters.A.size();
    long b_size = parameters.b.size();
    // print out all the parameters in (full, not sparse) matrix form
    // name of decision variables
    std::cout << "*******************************************************" << std::endl;;
    std::cout << "Decision Variables" << std::endl;
    for (int x_index = 0; x_index < x_size; ++x_index) {
        std::cout << "x" << x_index << " : " << parameters.x[x_index] << std::endl;
    }
    // print A
    std::cout << "*******************************************************" << std::endl;;
    std::cout << "A:" << std::endl;
    for (int row_index = 0; row_index < A_rowsize; ++row_index) {
        for (int x_index = 0; x_index < x_size; ++x_index) {
            // find if coefficient for current x is 0
            if (parameters.A[row_index].find(x_index) == parameters.A[row_index].end()) { // coefficient is 0
                std::cout << 0 << " ";
            }
            else {
                std::cout << parameters.A[row_index][x_index] << " ";
            }
        }
        std::cout << std::endl;
    }
    // print b
    std::cout << "*******************************************************" << std::endl;;
    std::cout << "b:" << std::endl;
    for (int row_index = 0; row_index < A_rowsize; ++row_index) {
        if (parameters.b.find(row_index) == parameters.b.end()) { // b[row_index] is 0
            std::cout << 0 << " ";
        }
        else {
            std::cout << parameters.b[row_index] << " ";
        }
        std::cout << std::endl;
    }
}

// test on functions in ioNonparametricModel
void inputModelTest() {
    std::string filePath = "/Users/sonny/Documents/predictorDB/model01.txt";
    twoStageParameters parameters = readTwoStageParameters(filePath);
    printTwoStageParameters(parameters);
}

void inputOneStageModelTest() {
    std::string filePath = "/Users/sonny/Documents/DebuggingPlayground/oneStagePiecewiseLP2/model.txt";
    oneStageParameters parameters = readOneStageParameters(filePath);
    printOneStageParameters(parameters);
}
