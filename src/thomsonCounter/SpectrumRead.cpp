#include "thomsonCounter/SpectrumRead.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

bool readSpectrumFromT(std::string filename, double &T0, double &dT, uint &N_T, darray &S_T, uint UN_USE_CHANNELS)
{
    S_T.clear();
    std::ifstream fin;
    fin.open(filename);

    if (fin.is_open())
    {
        std::string line;
        std::getline(fin, line);

        uint channels = 0;
        std::istringstream iss(line);
        std::string word;
        while (iss >> word)
            channels++;

        N_T = 0;

        while (!fin.fail())
        {
            double T;
            fin >> T;

            if (fin.fail())
                break;

            for (uint i = 0; i < channels; i++)
            {
                double value;
                fin >> value;
                S_T.push_back(value);
            }

            for (uint i = 0; i < UN_USE_CHANNELS; i++)
                S_T.push_back(0);

            if (N_T == 0)
            {
                T0 = T;
            }
            else if (N_T == 1)
            {
                dT = T - T0;
            }

            N_T++;            
        }
        
    }
    else
    {
        std::cerr << "не удалось открыть файл\n";
        return false;
    }

    fin.close();
    return true;
}