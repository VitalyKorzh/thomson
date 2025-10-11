#include "thomsonCounter/SRF.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>


bool readSRF(std::string filename, darray &SRF, double &lMin, double &lMax, double &dl, uint &N_LAMBDA, uint &N_CHANNELS, uint UN_USE_CHANNELS) 
{
    SRF.clear();

    std::ifstream fin;
    fin.open(filename);

    if (fin.is_open()) 
    {
        std::string line;
        std::getline(fin, line);

        std::istringstream iss(line);
        std::string word;
        uint channels = 0;
        while (iss >> word)
            channels++;

        channels-=2;
        N_CHANNELS = channels;

        N_LAMBDA = 0;

        double lambda, dontUse;
        darray SRF_T;

        while (!fin.fail())
        {
            fin >> lambda >> dontUse;

            if (fin.fail())
                break;


            for (uint i = 0; i < channels; i++)
            {
                double value;
                fin >> value;
                SRF_T.push_back(value);
            }

            for (uint i = 0; i < UN_USE_CHANNELS; i++)
                SRF_T.push_back(0);

            if (N_LAMBDA == 0)
                lMin = lambda;
            else if (N_LAMBDA == 1)
                dl = lambda - lMin;
            
            N_LAMBDA++;
        }

        lMax = lambda;
        N_CHANNELS+=UN_USE_CHANNELS;
        SRF.resize(SRF_T.size());

        for (uint i = 0; i < N_CHANNELS; i++)
        {
            for (uint j = 0; j < N_LAMBDA; j++)
                SRF[i*N_LAMBDA+j] = SRF_T[j*N_CHANNELS+i];
        }

    }
    else 
    {
        std::cerr << "не удалось открыть файл: " << filename << "\n";
        return false;
    }

    fin.close();
    return true;
}