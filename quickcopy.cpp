//
// Created by cornelia.blanke on 22.10.2024.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <direct.h>
using namespace std;


// https://stackoverflow.com/questions/9125122/how-to-copy-a-file-from-a-folder-to-another-folder
// copy in binary mode
void copyFile(const string& SRC, const string& DEST, bool check=true)
{
    std::ifstream src(SRC, std::ios::binary);
    if (src.is_open()) {
        std::ofstream dest(DEST, std::ios::binary);
        if (dest.is_open())
            dest << src.rdbuf();
    }
    else if (check) {   // print error message and exit(1)
        cerr << SRC << " not found" << endl;
        exit(1);
    }
}

int main()
{
    string OldDir, OldFiles, NewDir, NewFiles;
    //const string cases="C:\\Users\\cornelia.blanke\\OneDrive - HESSO\\Cornelia\\Enseignement\\24-25\\INTE\\Simulation\\Cases\\";
    const string cases="Cases\\";
    string line;

    // user entries
    cout << "Old Case Directory?" << endl;
    cin >> OldDir;
    cout << "Old Case Name?" << endl;
    cin >> OldFiles;
    cout << "New Case Directory? [" << OldDir << "]" << endl;
    cin.ignore();
    getline( cin, line);
    if ( !line.empty() ) {
        istringstream stream(line);
        stream >> NewDir;
        _mkdir((cases+NewDir).c_str());
    }
    else
        NewDir = OldDir;
    cout << "New Case Name?" << endl;
    cin >> NewFiles;

    const string src = cases+OldDir+"\\"+OldFiles;
    const string dest= cases+NewDir+"\\"+NewFiles;
    bool over=true;

    // what to do with existing file?
    if (FILE *file = fopen((dest+".csv").c_str(), "r")) {
        fclose(file);
        cout << "Override? (y/n) [n]" << endl;
        cin.ignore();
        getline( cin, line);
        if ( line[0]!='y' && line[0]!='Y' )
            over=false;
    }

    if (over){
        // copy files
        copyFile(src+".csv", dest+".csv");
        copyFile(src+"_diam.csv", dest+"_diam.csv");
        copyFile(src+"_length.csv", dest+"_length.csv");
        copyFile(src+"_pipes.csv", dest+"_pipes.csv");
        copyFile(src+"_SST.csv", dest+"_SST.csv");

        // copy files if exist
        copyFile(src+"_HP.csv", dest+"_HP.csv", false);
        copyFile(src+"_SST_dT.csv", dest+"_SST_dT.csv", false);
        copyFile(src+"_SST_demand.csv", dest+"_SST_demand.csv", false);

        // copy-modify settings file
        ifstream file(src+"_settings.txt", ios::in);
        ofstream outfile(dest+"_settings.txt");
        if (file.is_open() && outfile.is_open()) {
            while (getline(file, line)) {
                if (line.rfind("Case Directory", 0) == 0)
                    line.replace(line.find(OldDir), sizeof(OldDir) - 1, NewDir);
                else if (line.rfind("Case Name", 0) == 0)
                    line.replace(line.find(OldFiles), sizeof(OldFiles) - 1, NewFiles);

                outfile << line <<"\n";
            }
            file.close();
            outfile.close();
        }
    }
    return 0;
}
