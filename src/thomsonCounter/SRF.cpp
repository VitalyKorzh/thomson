#include "../../include/thomsonCounter/SRF.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>


bool readSRF(std::string filename, darray &SRF, double &lMin, double &lMax, double &dl, uint &N_LAMBDA, uint N_CHANNELS) 
{
    SRF.clear();
    lMin = 0.;
    lMax = 0.;
    dl = 0.;
    N_LAMBDA = 0;

    std::ifstream fin;
    fin.open(filename);

    if (fin.is_open()) 
    {

        std::string line;
        std::getline(fin, line);

        while (std::getline(fin, line))
        {
            std::istringstream iss(line);
            double l;
            iss >> l;

            if (iss.fail())
                break;
            
            double value;
            iss >> value;
            for (uint i = 0; i < N_CHANNELS; i++)
            {
                iss >> value;
                if (iss.fail())
                    value = 0.;

                SRF.push_back(value);
            }

            if (N_LAMBDA == 0)
                lMin = l;
            else if (N_LAMBDA == 1)
                dl = l - lMin;
            N_LAMBDA++;
            lMax = l;
        }
        
        darray SRF_TEMP(SRF.size());
        for (uint ch = 0; ch < N_CHANNELS; ch++)
        {
            for (uint it = 0; it  < N_LAMBDA; it++)
                SRF_TEMP[ch*N_LAMBDA+it] = SRF[ch + N_CHANNELS*it];
        }

        SRF = SRF_TEMP;

        // std::string line;
        // std::getline(fin, line);

        // std::istringstream iss(line);
        // std::string word;
        // uint channels = 0;
        // while (iss >> word)
        //     channels++;

        // channels-=2;
        // N_CHANNELS = channels;

        // N_LAMBDA = 0;

        // double lambda, dontUse;
        // darray SRF_T;

        // while (!fin.fail())
        // {
        //     fin >> lambda >> dontUse;

        //     if (fin.fail())
        //         break;


        //     for (uint i = 0; i < channels; i++)
        //     {
        //         double value;
        //         fin >> value;
        //         SRF_T.push_back(value);
        //     }

        //     for (uint i = 0; i < UN_USE_CHANNELS; i++)
        //         SRF_T.push_back(0);

        //     if (N_LAMBDA == 0)
        //         lMin = lambda;
        //     else if (N_LAMBDA == 1)
        //         dl = lambda - lMin;
            
        //     N_LAMBDA++;
        // }

        // lMax = lambda;
        // N_CHANNELS+=UN_USE_CHANNELS;
        // SRF.resize(SRF_T.size());

        // for (uint i = 0; i < N_CHANNELS; i++)
        // {
        //     for (uint j = 0; j < N_LAMBDA; j++)
        //         SRF[i*N_LAMBDA+j] = SRF_T[j*N_CHANNELS+i];
        // }

    }
    // else 
    // {
    //     std::cerr << "не удалось открыть файл: " << filename << "\n";
    //     return false;
    // }

    if (N_LAMBDA == 0)
    {
        N_LAMBDA = 2;
        lMin = 0.;
        lMax = 1064.;
        dl = lMax - lMin;
        SRF.resize(N_CHANNELS*2, 0.);
    }

    fin.close();
    return true;
}