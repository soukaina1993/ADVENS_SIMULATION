//
// Created by cornelia.blanke on 24.05.2023.
//

#include "network/QNetwork.h"
#include <chrono> //clock()

using namespace std;

int main(int argc, char *argv[]) {
    if (argc==2) {
        const char * inputfile = argv[1];
        clock_t start = clock();

        /* ------------------------------------------------------------------------------------------ */

        QNetwork Network(inputfile);
        Network.init();     // defines filepath, connects to database and downloads files
        Network.extract_quantum_network();

        /* ------------------------------------------------------------------------------------------ */

        //inputfile is now specified as command line argument:
        //QNetwork Network("C:/Users/cornelia.blanke/Documents/adv/Cases/Esta_NEW/Esta_NEW_settings.txt");
        //QNetwork Network(inputfile);
        Network.init_network();     // inits computation

        //while (Network.timestep < 1) {
        while (Network.timestep < Network.max_timestep) {
            Network.init_timestep();
            Network.solve_timestep();
            Network.save_results();
            Network.next_timestep();
        }

        Network.finalize(false);     // upload geopackages
        cout << endl << "finished" << endl;

        /* ------------------------------------------------------------------------------------------ */

        clock_t end = clock();
        cout << "Runtime: " << (double) (end - start) / CLOCKS_PER_SEC << " seconds" << endl;
    }
    else {
        cerr << "Please specify a settings file" << endl;
        exit(1);
    }
    return 0;
}