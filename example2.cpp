//
// Created by cornelia.blanke on 08.05.2024.
//

/*** before execution, Python must be added to the PATH variable, for example
 * set PATH=$PATH$;C:\Users\cornelia.blanke\AppData\Local\Programs\Python\Python310
 *
 * if packages need to be imported, add for example
 * set PYTHONPATH=C:\Users\cornelia.blanke\Documents\PyPipes\venv\Lib\site-packages
 * or call sys.path.append
 * ***/

#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <iostream>
#include <vector>

namespace py = pybind11;

int main() {
    //py::scoped_interpreter python;
    py::scoped_interpreter guard{};

    auto math = py::module::import("math");
    double root_two = math.attr("sqrt")(2.0).cast<double>();

    std::cout << "The square root of 2 is: " << root_two << "\n";

    // append pathes to modules
    py::exec(R"(
        import sys
        sys.path.append("C:\\Users\\cornelia.blanke\\Documents\\PyPipes\\venv\\Lib\\site-packages")
        sys.path.append("C:\\Users\\cornelia.blanke\\Documents\\PyPipes\\dhn")
        del sys
    )");

    // import numpy
    auto np = py::module::import("numpy");
    double e = np.attr("e").cast<double>();
    std::cout << "e is: " << e << "\n";
    py::array_t<double> rdata = np.attr("random").attr("rand")(10);
    double *carray = rdata.mutable_data();
    std::cout << "carray is: ";
    for (int i=0;i<10; i++){
        std::cout << carray[i] << " ";
    }
    std::cout << std::endl;

    //import own class
    //py::object pypipes = py::module::import("pypipes");

    // not working:
//    system(R"(
//        C:\Users\cornelia.blanke\Documents\PyPipes\venv\scripts\activate.bat; cd C:\Users\cornelia.blanke\Documents\PyPipes; python -m dhn.DHN_Sierre
//    )");


    return 0;
}