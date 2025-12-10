//
// Created by cornelia.blanke on 14.09.2023.
//
#include "externals/CoolProp/include/CoolProp.h"
#include <iostream>

int main()
{
    double p_V=4000000, T_1=450+273.15;
    std::cout << CoolProp::PropsSI("H","P",p_V,"T",T_1,"H2O") << std::endl;
    std::cout << CoolProp::PropsSI("S","P",p_V,"T",T_1,"H2O") << std::endl;
    return 1;
}