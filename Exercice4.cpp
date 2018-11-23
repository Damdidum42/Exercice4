#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "ConfigFile.tpp" // Fichier .tpp car inclut un template

using namespace std;

class Exercice4
{

private:
  double t, dt, tFin;
  double G,p0,RT;
  double mT,mL,mA;
  vector<double> y;
  double xT_0,yT_0,xL_0,yL_0,xA_0,yA_0,vTx_0,vTy_0,vLx_0,vLy_0,vAx_0,vAy_0;
  int sampling;
  int last;
  ofstream *outputFile;

  void printOut(bool force)
  {
    if((!force && last>=sampling) || (force && last!=1))
    {
      *outputFile << t << " "  << G << endl;
      last = 1;
    }
    else
    {
      last++;
    }
  }

  void step()
  {

  }


public:

  Exercice4(int argc, char* argv[])
  {
    string inputPath("configuration.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice3 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice3 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    tFin     = configFile.get<double>("tFin");
    G        = configFile.get<double>("G");
    p0       = configFile.get<double>("p0");
    RT       = configFile.get<double>("RT");

    mT       = configFile.get<double>("mT");
    mL       = configFile.get<double>("mL");
    mA       = configFile.get<double>("mA");

    xT_0     = configFile.get<double>("xT_0");
    yT_0     = configFile.get<double>("yT_0");
    vTx_0    = configFile.get<double>("vTx_0");
    vTy_0    = configFile.get<double>("vTy_0");

    xL_0     = configFile.get<double>("xL_0");
    yL_0     = configFile.get<double>("yL_0");
    vLx_0    = configFile.get<double>("vLx_0");
    vLy_0    = configFile.get<double>("vLy_0");

    xA_0     = configFile.get<double>("xA_0");
    yA_0     = configFile.get<double>("yA_0");
    vAx_0    = configFile.get<double>("vAx_0");
    vAy_0    = configFile.get<double>("vAy_0");

    dt       = configFile.get<double>("dt");
    sampling = configFile.get<int>("sampling");

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);
  };

  ~Exercice4()
  {
    outputFile->close();
    delete outputFile;
  };

  double acc(double t,double x, double v){
      //FONCTION D'ACCELERATION
  }

  void run()
  {
    t = 0.;
    last = 0;
    printOut(true);
    while( t < tFin-0.5*dt )
    {
      step();
      t += dt;
      printOut(false);
    }
    printOut(true);
  };

};


int main(int argc, char* argv[])
{
  Exercice4 engine(argc, argv);
  engine.run();
  return 0;
}
