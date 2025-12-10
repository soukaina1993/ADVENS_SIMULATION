//
// Created by cornelia.blanke on 28.08.2023.
//

/*
 * Inside MATLAB, install the MinGW-w64 Add-On as it needs to be the 64bit version!
 *
 * MATLAB function mex_test multiplies an array by 2
 * example: mex_test([1,2,3,4,5]) yields [2,4,6,8,10]
 */

#include "mex.hpp"
#include "mexAdapter.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;

class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {

        // Validate arguments
        checkArguments(outputs, inputs);

        // Implement function
        TypedArray<double> doubleArray = std::move(inputs[0]);
        for (auto& elem : doubleArray) {
            elem *= 2;
        }

        // Assign outputs
        outputs[0] = doubleArray;
    }

    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
        ArrayFactory factory;
        if (inputs[0].getType() != ArrayType::DOUBLE ||
            inputs[0].getType() == ArrayType::COMPLEX_DOUBLE)
        {
            matlabPtr->feval(u"error", 0,
                             std::vector<Array>({ factory.createScalar("Input must be double array") }));
        }

        if (outputs.size() > 1) {
            matlabPtr->feval(u"error", 0,
                             std::vector<Array>({ factory.createScalar("Only one output is returned") }));
        }
    }
};
