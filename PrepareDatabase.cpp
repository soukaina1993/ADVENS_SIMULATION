//
// Created by cornelia.blanke on 23.06.2025.
//

#include "network/QNetwork.h"
#include "db_connect/DBConnection.h"

int main(int argc, char *argv[]) {
    if (argc==2) {
        const char *inputfile = argv[1];
        QNetwork Network(inputfile);
        Network.set_filepath();
        Network.connect_to_database();
        Network.move_geopackage_to_database("Network", false);  // true: deletes local file
        Network.move_geopackage_to_database("Data", false);
    }
    else {
        cerr << "Please specify a settings file" << endl;
        exit(1);
    }
    return 0;
}