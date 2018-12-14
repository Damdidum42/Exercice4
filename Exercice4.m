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
dt= linspace(0.1,5,nsimul) % TODO: Choisir des valeurs de dt pour faire une etude de convergence

%EPSILON
epsilon = logspace(-7,-1,nsimul)




paramstr = 'dt'; % Nom du parametre a scanner (changer ici 'dt' ou 'Omega' ou autre)
param = dt; % Valeurs du parametre a scanner (changer ici dt ou Omega ou autre)

%% CONSTANTES DE SIMULATION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RT = 6378.1*1000 %rayon terrestre
RL = 1737 * 1000 %rayon lunaire
RA = 300 %rayon du vaisseau (pour qu'il soit visible)

h = 10000+RT
vmax_th = 11121.58

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

if strcmp(paramstr, 'dt')
    erreur_h = zeros(1,nsimul);
    erreur_v = zeros(1,nsimul);    
elseif strcmp(paramstr, 'epsilon')
    erreur_h = zeros(1,nsimul);
    erreur_v = zeros(1,nsimul);
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
    
    if strcmp(paramstr, 'dt') || strcmp(paramstr, 'epsilon')
        dist = ((xT-xA).^2 + (yT-yA).^2).^(0.5)
        erreur_h(i) = abs( min(dist) -h )
        v = (vAx.^2 + vAy.^2).^(0.5)
        erreur_v(i) = abs(max(v)-vmax_th)
    end

end

%% Figures %%
%%%%%%%%%%%%%

if strcmp(paramstr, 'dt')
    %figues pour l'étude de la convergence
    figure('name','alti convergence dt')
    loglog(dt,erreur_h,'+')
    grid on
    xlabel('log(\Delta t)')
    ylabel('distance minimale avec la surface terrestre [m]')
    
    figure('name','vitesse convergence dt')
    loglog(dt,erreur_v,'+')
    grid on
    xlabel('log(\Delta t)')
    ylabel('v_{max} [m/s]')
end

if strcmp(paramstr, 'epsilon')
    %figues pour l'étude de la convergence
    figure('name','alti convergence dt')
    loglog(epsilon,erreur_h,'+')
    grid on
    xlabel('log(\epsilon)')
    ylabel('distance minimale avec la surface terrestre [m]')
    
    figure('name','vitesse convergence dt')
    loglog(epsilon,erreur_v,'+')
    grid on
    xlabel('log(\epsilon)')
    ylabel('vitesse maximale [m/s]')
end