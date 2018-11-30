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

  double fact = 0.99; //pour éviter les boucles infinies
  int n = 4; //ordre de convergence du schéma utilisé

  valarray<double> y;
  valarray<double> mass;

  int sampling;
  int last;
  ofstream *outputFile;

  void printOut(bool force){ //ECRIT LES DONNEES DE LA SIMULATION DANS UN FICHIER DE SORTIE
    if((!force && last>=sampling) || (force && last!=1))
    {
        *outputFile << t;
        for(int i(0); i<y.size() ; i++){
            *outputFile << " " << y[i];
        }
      //*outputFile << t << " " << xT_0 << " " << yT_0 << " " << xL_0 << " " << yL_0 << " " << xA_0 << " " << yA_0 << " " << vTx_0 << " " << vTy_0 << " " << vLx_0 << " " << vLy_0 << " " << vAx_0 << " " << vAy_0 << endl;

      //ECRIRE LES ENERGIES
      *outputFile<< " " << Emec(y) << Ptot(y)[0] << Ptot(y)[1];

      *outputFile << endl;

      last = 1;
    }
    else
    {
      last++;
    }
  }

  void step(){//EFFECTUE UN PAS DE TEMPS DE LA SIMULATION
    dt_adapt(y);
    y = do_one_step(y,dt);
  }

  valarray<double> do_one_step(valarray<double> y, double dt){ //SIMULE UN PAS DE TEMPS DE LA SIMULATION (POUR LE CALCUL DU PAS DE TEMPS
    //schéma de Runge Kutta d'ordre 4:
    valarray<double> k1 (y.size());
    valarray<double> k2 (y.size());
    valarray<double> k3 (y.size());
    valarray<double> k4 (y.size());

    k1 = dt*f(t,y);
    k2 = dt*f(t+0.5*dt,y+0.5*k1);
    k3 = dt*f(t+0.5*dt,y+0.5*k2);
    k4 = dt*f(t+dt,y+k3);

    return y + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4);
  }

  valarray<double> f(double t, valarray<double> y){ //FONCTION D'ACCELERATION DU VECTEUR Y

    valarray<double> acc(y.size());
    size_t nbcorps(mass.size());
    valarray<double> distanceij(2); // distance entre 2 vecteurs position
    double norme(0); // norme d'une difference de vecteur position
    double norme3(0);
    // dy/dt = v

    acc[slice(0,2*nbcorps,1)]= y[slice(2*nbcorps,2*nbcorps,1)];

    //dv/dt = a
	  for (size_t i(0); i<nbcorps;i++){
        for(size_t j(0);j<nbcorps;j++){

            if(i!=j){
                valarray<double> ri(y[slice(i*2,2,1)]);
                valarray<double> rj(y[slice(j*2,2,1)]);
                distanceij=ri-rj;

                norme=norm(distanceij);
                norme3= pow(norme,3);
				acc[slice(2*i+2*nbcorps,2,1)]+= G*mass[j]*distanceij/norme3;
            }
        }
	  }
        return acc;

  }

  double rho(double r){//FONCTION QUI RETOURNE LA DENSITE EN FONCTION DE LA DISTANCE r
    return p0*exp(-(r-RT)/lambda);
  }

  double Emec(valarray<double> y){//CALCUL DE l'ENERGIE MECANIQUE DU SYSTEME
    double Ecin(0.);
    double Epot(0.);

    valarray<double> posT = y[slice(0,2,1)];
    valarray<double> posL = y[slice(2,2,1)];
    valarray<double> posA = y[slice(4,2,1)];

    valarray<double> vT = y[slice(6,2,1)];
    valarray<double> vL = y[slice(8,2,1)];
    valarray<double> vA = y[slice(10,2,1)];

        Ecin = mass[0] * pow(norm(vT),2) + mass[1] * pow(norm(vL),2) + mass[2] * pow(norm(vA),2);
        Ecin *= 0.5;

        Epot = (mass[0]*mass[1])/norm(posT-posL) + (mass[0]*mass[2])/norm(posT-posA) + (mass[1]*mass[2])/norm(posL-posA);
        Epot *= -G;

    return Ecin + Epot;

  }

  valarray<double> Ptot(valarray<double> y){//CALCUL DE LA QUANTITE DE MOUVEMENT TOTALE DU SYSTEME
    valarray<double> P;
    valarray<double> vT = y[slice(6,2,1)];
    valarray<double> vL = y[slice(8,2,1)];
    valarray<double> vA = y[slice(10,2,1)];

    P = mass[0]*vT + mass[1]*vL + mass[2]*vA;
  }

  void dt_adapt(valarray<double> y){ //FONCTION POUR ADAPTER LE PAS DE TEMPS
    double new_dt;
    valarray<double> y1 = do_one_step(y,dt);
    valarray<double> y2 = do_one_step(do_one_step(y,dt*0.5),0.5*dt);
    double d  = norm(y1-y2);

    if(d > epsilon){
        while(d > epsilon){
            new_dt = dt*fact*pow( (epsilon/d) , 1.0/(n+1) );
            y1 = do_one_step(y,new_dt);
            y2 = do_one_step(do_one_step(y,new_dt*0.5),0.5*new_dt);
            d  = norm(y1-y2);
        }
        dt= new_dt;
    }

    else{
        dt = dt*fact*pow(( epsilon/d), 1.0/(n+1) );
    }
  }

public:

  Exercice4(int argc, char* argv[]){
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

    mass.resize(3);

    mass[0]       = configFile.get<double>("mT");
    mass[1]       = configFile.get<double>("mL");
    mass[2]       = configFile.get<double>("mA");

    y.resize(12);

    y[0] = configFile.get<double>("xT_0");
    y[1] = configFile.get<double>("yT_0");
    y[2] = configFile.get<double>("xL_0");
    y[3] = configFile.get<double>("yL_0");
    y[4] = configFile.get<double>("xA_0");
    y[5] = configFile.get<double>("yA_0");

    y[6] = configFile.get<double>("vTx_0");
    y[7] = configFile.get<double>("vTy_0");
    y[8] = configFile.get<double>("vLx_0");
    y[9] = configFile.get<double>("vLy_0");
    y[10]= configFile.get<double>("vAx_0");
    y[11]= configFile.get<double>("vAy_0");

    dt       = configFile.get<double>("dt");
    sampling = configFile.get<int>("sampling");

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);
  };

  ~Exercice4(){
    outputFile->close();
    delete outputFile;
  };

  void run(){
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

double norm(valarray<double> x){//RETOURNE LA NORME D'UN VECTEUR
double sum(0);
    for(int i(0); i<=x.size(); i++)
        sum += x[i]*x[i];
return sqrt(sum);
}
