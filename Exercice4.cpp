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
  double G,p0,RT,RL,lambda,C_x,d;

  double fact = 0.999; //pour éviter les boucles infinies
  int n = 4; //ordre de convergence du schéma utilisé

  valarray<double> y;
  valarray<double> mass;

  int sampling;
  int last;
  ofstream *outputFile;

  void printOut(){ //ECRIT LES DONNEES DE LA SIMULATION DANS UN FICHIER DE SORTIE
        *outputFile << t;
        for(int i(0); i<y.size() ; i++){
            *outputFile << " " << y[i];
        }
        *outputFile << endl;
  }

  void step(){//EFFECTUE UN PAS DE TEMPS DE LA SIMULATION
    cout<<"t= "<<t<<" ; y = (";
    for(int i(0); i<y.size(); i++)
    cout<<y[i]<<" , ";
    cout<<endl;
    dt_adapt(y);
    y = do_one_step(y,dt);

    valarray<double> pos_Appollo = y[slice(4,2,1)];
    crash((pos_Appollo));
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
				acc[slice(2*i+2*nbcorps,2,1)]+= -G*mass[j]*distanceij/norme3;
            }
        }
	  }


      /*cout<<"acc = (";
      for(int i(0); i<acc.size(); i++)
        cout<<acc[i]<<" , ";
      cout<<endl;*/


	  //ajout de la force de trainée au vaisseau Appollo
	  double S = M_PI*pow(d/2,2);
	  valarray<double> pos_terre = y[slice(0,2,1)];
	  valarray<double> pos_Appollo = y[slice(4,2,1)];
	  valarray<double> vit_terre = y[slice(6,2,1)];
	  valarray<double> vit_Appollo = y[slice(10,2,1)];

       //cout<<"terre = ("<<pos_terre[0]<<" , "<<pos_terre[1]<<")"<<" Appollo = ("<<pos_Appollo[0]<<" , "<<pos_Appollo[1]<<")"<<endl;
       //cout<<"terre-Appollo = "<<norm(pos_terre-pos_Appollo)<<endl;
	  valarray<double> F_train = -0.5*rho(norm(pos_terre-pos_Appollo))*S*C_x*norm(vit_terre-vit_Appollo)*(vit_Appollo-vit_terre);
      acc[slice(10,2,1)] += F_train;

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
    double new_dt(dt);
    valarray<double> y1 = do_one_step(y,dt);
    valarray<double> y_prime = do_one_step(y,dt*0.5);
    valarray<double> y2 = do_one_step(y_prime,0.5*dt);
    double d  = norm(y1-y2);

    cout<<"ADAPATION!"<<endl;
    cout<<"y1 = (";
    for(int i(0); i<y1.size(); i++)
    cout<<y1[i]<<" , ";
    cout<<endl;
    cout<<"y2 = (";
    for(int i(0); i<y2.size(); i++)
    cout<<y2[i]<<" , ";
    cout<<endl;
    cout<<"d= "<<d<<" ; epsilon = "<<epsilon<<endl;

    if(d >=  epsilon){
        while(d > epsilon){
            new_dt = dt*fact*pow( (epsilon/d) , 1.0/(n+1) );
            y1 = do_one_step(y,new_dt);
            y2 = do_one_step(do_one_step(y,new_dt*0.5),0.5*new_dt);
            d  = norm(y2-y1);
            cout<<"dt: "<<dt<<" -> "<< new_dt <<endl;
            dt= new_dt;
        }

    }

    else{
        dt = pow(( epsilon/d), 1.0/(n+1) );
        cout<<"dt reste:"<<dt <<endl;
    }
  }

  valarray<double> C_masse(valarray<double> masses, valarray<double> y){//calcule le centre de masse du systême terre lune
      valarray<double> C_m;
      double masse_tot(0);
      C_m.resize(2);
      for(int i(0); i<masses.size(); i++){
        masse_tot += masses[i];
      }

      for(int i(0); i<masses.size(); i++){
        valarray<double> element = y[slice(2*i,2,1)];
        C_m += masses[i]*element;
      }

      cout<<"C_m = ("<<C_m[0]<<" "<<C_m[1]<<")"<<endl;
      return C_m*(1/masse_tot);
  }

  valarray<double> V_masse(valarray<double> masses, valarray<double> y){//calcule le centre de masse du systême terre lune
      valarray<double> V_m;
      double masse_tot(0);
      size_t nbcorps = masses.size();
      V_m.resize(2);
      for(int i(0); i<nbcorps; i++){
        masse_tot += masses[i];
      }

      for(int i(0); i<nbcorps; i++){
        valarray<double> element = y[slice( 2*nbcorps +2*i,2,1)];
        V_m += masses[i]*element;
      }

      cout<<"V_m = ("<<V_m[0]<<" "<<V_m[1]<<")"<<endl;
      return V_m*(1/masse_tot);
  }

  void crash(valarray<double> r1){ //fonction qui nous permet de savoir si le vaisseau Appollo est entré en collision
    valarray<double> posTerre = y[slice(0,2,1)];
    valarray<double> posLune = y[slice(2,2,1)];

    if(norm(r1-posTerre)<= RT){
        cout<<"crash avec la terre"<<endl;
        exit(0);
    }
    if(norm(r1-posLune)<= RL){
      cout<<"crash avec la terre"<<endl;
      exit(0);
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
    tFin *= 60*60*24;//conversion jour->secondes
    epsilon  = configFile.get<double>("epsilon");
    G        = configFile.get<double>("G");
    p0       = configFile.get<double>("p0");
    RT       = configFile.get<double>("RT");
    RL       = configFile.get<double>("RL");
    lambda   = configFile.get<double>("lambda");
    C_x      = configFile.get<double>("C_x");
    d        = configFile.get<double>("d");

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

    //translation vers le centre de masse
    valarray<double> y_G = C_masse(mass,y);
    valarray<double> v_G = V_masse(mass,y);

    //passe au référentiel de centre de masse
    for(int i(0); i< mass.size(); i++){
        y[slice(2*i,2,1)] -=  y_G;
    }

    for(int i(0); i< mass.size(); i++){
        y[slice(2*mass.size() + 2*i,2,1)] -=  v_G;
    }

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);

       cout<<"ADAPATION!"<<endl;
    cout<<"y_initiale = (";
    for(int i(0); i<y.size(); i++)
    cout<<y[i]<<" , ";
    cout<<endl;

  };

  ~Exercice4(){
    outputFile->close();
    delete outputFile;
  };

  void run(){
    t = 0.0;
    printOut();
    while(t < tFin-(0.5*dt) )
    {
      //cout<<"t= "<<t<<" tFin = "<<tFin<<endl;
      step();
      t += dt;
      printOut();
    }
    printOut();
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
