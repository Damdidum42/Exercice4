
repertoire = 'C:\Users\damie\Desktop\EPFL-PH-BA2\PH-num\exercices\4\'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice4'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

%loading
filename = 'out/a.out'
data = load(filename)

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
%==============================================================
%constantes de la simulation

RT = 6378.1*1000 %rayon terrestre
RL = 1737 * 1000 %rayon lunaire
RA = 300 %rayon du vaisseau (pour qu'il soit visible)


%===============================================================
%figures

    figure('name',filename)
    plot(xT,yT)
    circle(xT(end),yT(end),RT)
    hold on
    %circle(xL(end),yL(end),RL)
    %plot(xL,yL)
    hold on 
    plot(xA,yA)
    circle(xA(end),yA(end),RT)
    axis equal
    xlabel('x [m]')
    ylabel('y [m]')
    
function h = circle(x,y,r)
hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
plot(xunit, yunit);
hold off
end