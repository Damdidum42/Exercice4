% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
% 
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.
%

%% Parametres %%
%%%%%%%%%%%%%%%%

%constantes de la simulation

repertoire = 'C:\Users\damie\Desktop\EPFL-PH-BA2\PH-num\exercices\4\'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice4'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

nsimul = 20; % Nombre de simulations a faire

%NSTEPS
dt= linspace(1,20,nsimul) % TODO: Choisir des valeurs de dt pour faire une etude de convergence

%EPSILON
epsilon = logspace(-7,-1,nsimul)




paramstr = 'epsilon'; % Nom du parametre a scanner (changer ici 'dt' ou 'Omega' ou autre)
param = epsilon; % Valeurs du parametre a scanner (changer ici dt ou Omega ou autre)

%% CONSTANTES DE SIMULATION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RT = 6378.1*1000 %rayon terrestre
RL = 1737 * 1000 %rayon lunaire
RA = 300 %rayon du vaisseau (pour qu'il soit visible)

G = 6,67408*10^(-11) 
mT =  5.972*10^(24)

h = 10000
r_0 =  314159000
v_0 = 1200
vmax_th = 11121.5838699

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = ['out/',paramstr, '=', num2str(param(i)), '.out'];
    %Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

if strcmp(paramstr, 'dt') || strcmp(paramstr, 'epsilon')
    erreur_h = zeros(1,nsimul);
    erreur_v = zeros(1,nsimul); 
    pas_de_temps = zeros(1,nsimul);
    Pmax = zeros(1,nsimul);
    Accmax = zeros(1,nsimul);
end

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    
    t= data(:,1);
    
    xT = data(:,2);
    yT = data(:,3);
    xL = data(:,4);
    yL = data(:,5);
    xA = data(:,6);
    yA = data(:,7);
    
    vTx = data(:,8);
    vTy = data(:,9);
    vLx = data(:,10);
    vLy = data(:,11);
    vAx = data(:,12);
    vAy = data(:,13);
    
    dt = data(:,14);
    AccApp = data(:,15);
    nbr_pas = data(:,16);
    Pfrott = data(:,17);
    Pfrott(1)=0
      
    if strcmp(paramstr, 'dt') || strcmp(paramstr, 'epsilon')
        
        %interpolation parabolique (essai)
        dist = ((xT-xA).^2 + (yT-yA).^2).^(0.5)
        [P,k,pk]=peak(t,dist)
        v = (vAx.^2 + vAy.^2).^(0.5)
        [V,j,vk] = peak(t,-v)
 
        %calculs
        erreur_h(i) = abs(pk-h-RT)
        erreur_v(i) = abs(-vk-vmax_th)
        pas_de_temps(i) = nbr_pas(end)
        
        Accmax(i) = max(AccApp)
        figure('name','Acc')
        plot(t,AccApp)
        Pmax(i) = min(Pfrott)
        figure('name','Pfrott')
        plot(t,Pfrott)
    end
    
    if strcmp(paramstr, 'epsilon')
        
        
    end

end

%% Figures %%
%%%%%%%%%%%%%

if strcmp(paramstr, 'dt')
    %figues pour l'étude de la convergence
%     figure('name','alti convergence dt')
%     fitloglog(dt,erreur_h)
%     grid on
%     xlabel('log(\Delta t)','FontSize', 14)
%     ylabel('log(|h_{min}-h|)','FontSize', 14)
    
    figure('name','alti convergence nbr_pas')
    fitloglog(pas_de_temps,erreur_h)
    grid on
    xlabel('log(#pas de temps)','FontSize', 14)
    ylabel('log(|h_{min}-h|)','FontSize', 14)
    
%     figure('name','vitesse convergence dt')
%     fitloglog(dt,erreur_v)
%     grid on
%     xlabel('log(\Delta t)','FontSize', 14)
%     ylabel('log(|v_{max}-v_{th}|)','FontSize', 14)
end

if strcmp(paramstr, 'epsilon')
    %figues pour l'étude de la convergence
%     figure('name','alti convergence epsilon')
%     fitloglog(epsilon,erreur_h)
%     grid on
%     xlabel('log(\epsilon)','FontSize', 14)
%     ylabel('log(|h_{min}-h|)','FontSize', 14)
    
    figure('name','alti convergence nbr_pas')
    fitloglog(pas_de_temps,erreur_h)
    grid on
    xlabel('log(#pas de temps)','FontSize', 14)
    ylabel('log(|h_{min}-h|)','FontSize', 14)
    
    figure('name','Pmax epsilon')
    loglog(epsilon,Pmax)
    grid on
    xlabel('log(\epsilon)','FontSize', 14)
    ylabel('log(P_{max})','FontSize', 14)
    
%     figure('name','vitesse convergence epsilon')
%     fitloglog(epsilon,erreur_v)
%     grid on
%     xlabel('log(\epsilon)')
%     ylabel('log(|v_{max}-v_{th}|)','FontSize', 14)
end

function [P,i,pik] = peak(X,Y)
[a,i]=min(Y)
p=polyfit(X(i-2:i+2),Y(i-2:i+2),2)
L = linspace(X(i-2),X(i+2),100)
P=polyval(p,L)
pik= min(P)
end

function fitloglog(x,y)
X=log10(x); Y=log10(y);  % convert both variables to log's
b=polyfit(X,Y,1);        % estimate coefficients
yhat=10.^[polyval(b,[X(1) X(end)])];  % evaluate at end points
loglog(x,y,'+')
hold on
loglog([x(1) x(end)],yhat) % add fitted line
theString = sprintf('y = %.3f x + %.3f', b(1), b(2));
text(1,1, theString, 'FontSize', 18,'unit','normalized');
end
