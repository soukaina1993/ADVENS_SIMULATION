//
// Created by lucile.schulthe on 01.10.2021.
//

#include "iofilecsv.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>

using namespace std;




void create(char* file_input)
{
    // file pointer
    fstream fout;

    // opens an existing csv file or creates a new file.
    fout.open(file_input, ios::out | ios::app);

    cout << "Enter the details of 5 students:"
         << " roll name maths phy chem bio"<< endl;

    int i, roll, phy, chem, math, bio;
    string name;

    // Read the input
    for (i = 0; i < 5; i++) {

        cin >> roll
            >> name
            >> math
            >> phy
            >> chem
            >> bio;

        // Insert the data to file
        fout << roll << ", "
             << name << ", "
             << math << ", "
             << phy << ", "
             << chem << ", "
             << bio
             << "\n";
    }
}

void read_record2(char* file_input)
{
    // File pointer
    cout<<"on est ici"<<endl;
    fstream fin;

    // Open an existing file
    fin.open(file_input, ios::in);

    // Get the roll number
    // of which the data is required
    int count = 0;


    // Read the Data from the file
    // as String Vector
    vector<string> row;
    string line, word, temp;

    while (fin >> temp) {
        row.clear();

        getline(fin, line);

        stringstream s(line);

        while (getline(s, word, ',')) {
            // add all the column data
            // of a row to a vector
            row.push_back(word);
        }

    }
    if (count == 0)
        cout << "Record not found\n";
}

