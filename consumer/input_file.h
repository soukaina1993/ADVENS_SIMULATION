//
// Created by lucile.schulthe on 28.01.2022.
//

#ifndef MAIN_CPP_INPUT_FILE_H
#define MAIN_CPP_INPUT_FILE_H

#include "sst.h"

class input_file: public sst {
public:
    input_file(vector <string> &data_input);
    void calcul_power(environment &ext){};

    void calcul_flow(double Tinlet) { };

    void init(){
        energy=0;
    }


};


#endif //MAIN_CPP_INPUT_FILE_H
