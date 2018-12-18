
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
    
    dt = data(:,14);
    AccApp = data(:,15);
    nbr_pas = data(:,16);
    Pfrott = data(:,17);
    Pfrott(1)=0
%==============================================================
%constantes de la simulation

RT = 6378.1*1000 %rayon terrestre
RL = 1737 * 1000 %rayon lunaire
RA = 30000 %rayon du vaisseau (pour qu'il soit visible)


%===============================================================
%figures

    figure('name',filename)
    plot(xT,yT)
    hold on
    %circle(xT(end),yT(end),RT)
    hold on
    circle(xT(end),yT(end),RT+10000)
%     hold on
%     plot(xL,yL)
%     circle(xL(end),yL(end),RL)
    hold on 
    plot(xA,yA)
    %circle(xA(end),yA(end),RT)
    axis equal
    xlabel('x [m]')
    ylabel('y [m]')

figure('name','bla')
r = ( (xA-xT).^2 + (yA-yT).^2 ).^(0.5)
%[P,i,pk]=peak(t,r)
plot(t,r)
xlabel('t [s]','FontSize', 14)
ylabel('r [m]','FontSize', 14)
% hold on
% L = linspace(t(i-2),t(i+2),100)
% plot( L ,P)

% figure('name','dt')
% plot(t,dt)
% xlabel('t [s]','FontSize', 14)
% ylabel('dt [s]','FontSize', 14)
% 
% figure('name','dt/r')
% plot(r,dt)
% xlabel('r [m]','FontSize', 14)
% ylabel('dt [s]','FontSize', 14)

figure('name','norme acc')
plot(t,AccApp)
xlabel('t')
ylabel('||acc||')

figure('name','Pfrott')
plot(t,Pfrott)
xlabel('t')
ylabel('P')

    
function h = circle(x,y,r)
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
plot(xunit, yunit);
end

function [P,i,pik] = peak(X,Y)
[a,i]=min(Y)
p=polyfit(X(i-2:i+2),Y(i-2:i+2),2)
L = linspace(X(i-2),X(i+2),100)
P=polyval(p,L)
pik= min(P)
end
