#include <iostream>
#include <vector>
#include <fstream>
#include <string>


uint findPoint(const std::vector<double> &t, double t0)
{
    uint k = 0;
    for (uint i = 0; i < t.size(); i++)
    {
        if (t[i] > t0)
        {
            k = i;
            break;
        }
    }
    return k;
}

void format(const char *file_time, const char *file_ranges_time, const char *file_ranges_points)
{
    const uint timeSize = 1000;
    std::vector <double> t(timeSize);

    std::ifstream fin;

    fin.open(file_time);

    for (uint i = 0; i < timeSize; i++)
        fin >> t[i];

    fin.close();


    const uint N_SP = 6;
    const uint N_CH = 8;


    std::ofstream fout;
    fout.open(file_ranges_points);

    fin.open(file_ranges_time);


    std::vector <double> t1(N_CH*N_SP);
    std::vector <double> t2(N_CH*N_SP);
    std::string line;
    std::getline(fin, line);
    for (uint ch = 0; ch < N_CH; ch++)
    {
        for (uint sp = 0; sp < N_SP; sp++)
        {
            fin >> t1[sp*N_CH+ch] >> t2[sp*N_CH+ch];
            t2[sp*N_CH+ch] += t1[sp*N_CH+ch];
        }
    }


    for (uint sp = 0; sp < N_SP; sp++)
    {
        fout << "sp" << sp << "\n";
        for (uint ch = 0; ch < N_CH; ch++)
        {
            fout << "\tch" << ch << "\n";
            if (t2[sp*N_CH+ch] > t1[sp*N_CH+ch]) {
                fout << "\t\t" << 0 << " " << findPoint(t, t1[sp*N_CH+ch]) << "\n";
                fout << "\t\t" << 0 << " " << 0 << "\n";
                fout  << "\t\t"<< findPoint(t, t2[sp*N_CH+ch]) << " " << 1 << "\n";
                fout  << "\t\t" << findPoint(t, t1[sp*N_CH+ch]) << "\n";
            }
            else {
                fout << "\t\t" << 0 << " " << findPoint(t, t1[sp*N_CH+ch]) << "\n";
                fout << "\t\t" << 0 << " " << 0 << "\n";
                fout  << "\t\t"<< 999 << " " << 1 << "\n";
                fout  << "\t\t" << findPoint(t, t1[sp*N_CH+ch]) << "\n";
            }
            fout << "\t\t0.02 10 10\n";
        }
    }



    fin.close();
    fout.close();
}