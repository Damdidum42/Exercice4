#include <iostream>
#include <fstream>
#include <cmath>
#include <valarray>
#include <iomanip>
#include "ConfigFile.tpp" // Fichier .tpp car inclut un template

using namespace std;

//--------------------METHODE POUR GERER LES VALARRAY----------
double norm(valarray<double>x);

class Corps{
private:
    double masse;
    double x_0;
    double y_0;
    double vx_0;
    double vy_0;

public:
    Corps(double m, double x, double y, double vx, double vy)
    :masse(m),x_0(x),y_0(y),vx_0(vx),vy_0(vy)
    {}

};

class Exercice4
{

private:
  double t, dt, tFin, epsilon;
  double G,p0,RT,lambda;
  double mT,mL,mA;

  valarray<double> y;

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
    //mettre un schéma de Runge Kutta d'ordre 4 ici
    valarray<double> k1 (y.size());
    valarray<double> k2 (y.size());
    valarray<double> k3 (y.size());
    valarray<double> k4 (y.size());
  }

  valarray<double> do_one_step(valarray<double> y, double dt){
    //faire un seul pas de temps
  }

  valarray<double> f(double t, valarray<double> y){ //FONCTION D'ACCELERATION
    valarray<double> f (y.size());
    /*
    */

    f[0] = y[2]; // d/dt xT = vTx
    f[1] = y[3] ;// d/dt yT = vTy
  }

  valarray<double> F_grav_ij(double mi, valarray<double> ri, double mj, valarray<double> rj)
  {
  }
  double rho(double r){
    return p0*exp(-(r-RT)/lambda);
  }

  double adapt(valarray<double> y1, valarray<double> y2, double Delta_t){

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
    epsilon  = configFile.get<double>("epsilon");
    G        = configFile.get<double>("G");
    p0       = configFile.get<double>("p0");
    RT       = configFile.get<double>("RT");
    lambda   = configFile.get<double>("lambda");

    mT       = configFile.get<double>("mT");
    mL       = configFile.get<double>("mL");
    mA       = configFile.get<double>("mA");

    y.resize(12);

    y[0] = configFile.get<double>("xT_0");
    y[1] = configFile.get<double>("yT_0");
    y[2] = configFile.get<double>("vTx_0");
    y[3] = configFile.get<double>("vTy_0");

    y[4] = configFile.get<double>("xL_0");
    y[5] = configFile.get<double>("yL_0");
    y[6] = configFile.get<double>("vLx_0");
    y[7] = configFile.get<double>("vLy_0");

    y[8] = configFile.get<double>("xA_0");
    y[9] = configFile.get<double>("yA_0");
    y[10]= configFile.get<double>("vAx_0");
    y[11]= configFile.get<double>("vAy_0");

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

//--------------------------FONCTION POUR GÈRER LES VECTEURS-------------------------

double norm(valarray<double> x){
double sum(0);
    for(int i(0); i<=x.size(); i++)
        sum += x[i]*x[i];
return sqrt(sum);
}
