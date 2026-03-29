#include "thomsonCounter/SpectrumRead.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#define DEFAULT_T0 100.

bool readSpectrumFromT(std::string filename, double &T0, double &dT, uint &N_T, darray &S_T, uint N_CHANNELS)
{
    std::ifstream fin;
    fin.open(filename);

    S_T.clear();
    N_T = 0;
    T0 = 0.;
    dT = 0.;

    if (fin.is_open())
    {
        std::string line;
        std::getline(fin, line);

        while (std::getline(fin, line))
        {
            std::istringstream iss(line);
            double T;

            iss >> T;

            if (iss.fail())
                break;

            for (uint i = 0; i < N_CHANNELS; i++)
            {
                double value;
                iss >> value;
                if (iss.fail())
                    value = 0;

                S_T.push_back(value);
            }

            if (N_T == 0)
                T0 = T;
            else if (N_T == 1)
                dT = T - T0;
            N_T++;
        }
        
        darray S_T_TEMP(S_T.size());
        for (uint ch = 0; ch < N_CHANNELS; ch++)
        {
            for (uint it = 0; it < N_T; it++)
                S_T_TEMP[ch*N_T+it] = S_T[ch + N_CHANNELS*it];
        }
        S_T = S_T_TEMP;

        //uint channels = 0;
        // std::istringstream iss(line);
        // std::string word;
        // while (iss >> word)
        //     channels++;

        // N_T = 0;

        // while (!fin.fail())
        // {
        //     double T;
        //     fin >> T;

        //     if (fin.fail())
        //         break;

        //     for (uint i = 0; i < channels; i++)
        //     {
        //         double value;
        //         fin >> value;
        //         S_T.push_back(value);
        //     }

        //     // for (uint i = 0; i < UN_USE_CHANNELS; i++)
        //     //     S_T.push_back(0);

        //     if (N_T == 0)
        //     {
        //         T0 = T;
        //     }
        //     else if (N_T == 1)
        //     {
        //         dT = T - T0;
        //     }

        //     N_T++;            
        // }


        // darray S_T_TEMP(S_T.size());

        // for (uint ch = 0; ch < channels+UN_USE_CHANNELS; ch++)
        // {
        //     for (uint it = 0; it < N_T; it++)
        //         S_T_TEMP[ch*N_T+it] = S_T[ch+(channels+UN_USE_CHANNELS)*it];
        // }
        // S_T = S_T_TEMP;
        
    }
    // else
    // {
    //     std::cerr << "не удалось открыть файл: "  << filename << "\n";
    //     return false;
    // }

    if (N_T == 0) //если нет данных значит начальная темпераутра 100 эВ для всех
    {
        N_T = 1;
        T0 = DEFAULT_T0;
        dT = 0.;
        S_T.resize(N_CHANNELS, 0.);
        for (uint ch = 0; ch < N_CHANNELS; ch++)
            S_T[ch] = 1.;
    }

    fin.close();
    return true;
}