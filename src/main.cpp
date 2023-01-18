#include "sasa.hpp"

int main(int argc, char* argv[])
{
    srand(time(NULL));
    SASA *sasa = new SASA(argv[1]);
    sasa->run();
    delete sasa;
    return 0;
}
