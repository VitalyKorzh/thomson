#include <iostream>
#include "ThomsonGUI.h"


int main(int argc, char**argv) 
{
    ThomsonGUI thomsonGUI(gClient->GetRoot(), 600, 550, new TApplication("app", &argc, argv));
    thomsonGUI.run();

    return 0;
}