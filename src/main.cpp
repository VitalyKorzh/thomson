#include <iostream>
#include "ThomsonGUI.h"


int main(int argc, char**argv) 
{
    ThomsonGUI thomsonGUI(gClient->GetRoot(), 550, 600, new TApplication("app", &argc, argv), "Thomson", "thomson", 1064., 1000, 48, 11, 6, 8, 2, 7, 3, 6, 1, 10000);
    thomsonGUI.run();

    return 0;
}