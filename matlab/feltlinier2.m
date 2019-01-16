function feltlinier2
%% Dette program beregner det elektriske felt fra en spredende partikel.
clear all; clc; close all
%Definitioner
global e1 e2 b1 b2 c m Beta gamma t0 d
e1=1.6*2.99*1E-10; %elektron
e2=1.6*2.99*25E-10; %Spredningspartikel
b1=0.2E-10; %Impact parameter på første partikel
b2=.5E-10; %Impact parameter på anden partikel
c=2.99E10; % Lysets hastighed i cm/s
m=9.11E-28; %Elektronens masse
d=20E-10; %Afstand mellem stødpartikler

T0=0.5E-18;

Beta=0.9; %Beta faktoren
gamma=sqrt(1-Beta^2)^(-1); %Gamma faktoren

t0=1E-18; %Tidspunktet for billedet af feltlinierne
figure
%% Løs af differentialligning for forskellige startbetingelser

sspan=[0 1E-7]; %Løsningsområde for difflign
options=odeset('Abstol',1E-15,'Reltol',3E-14); %Indstillinger for ode45
nframes = 200;
n=20; % Det halve antal feltlinier
hold on
set(gca,'PlotBoxAspectRatio',[1 1 1])

for ti=1:nframes
    
    clf
    hold on
    set(gca,'PlotBoxAspectRatio',[1 1 1])
    
    t0 = 0.5E-18 * (1.5 * ti - 0.5 * nframes)/nframes; %Tidspunktet for billedet af feltlinierne
    
    for N = 0:n-1
        
        r0 = pos2(t0) + 0.5 * b1 * [cos(atan(gamma * tan(pi * N / n))); sin(atan(gamma*tan(pi*N/n))); 0]'; %Startbetingelse
        [s,r]=ode23(@difflign,sspan,r0,options); %Løsning af differentialligning
        
        plot(r(:,1),r(:,2),'b') %Plotning af løsning

        if N~=n/2
            r0=pos2(t0)+.5*b1*[-cos(atan(gamma*tan(pi*N/n))); sin(atan(gamma*tan(pi*N/n))); 0]'; %Startbetingelse
            [s,r]=ode23(@difflign,sspan,r0,options); %Løsning af differentialligning
            
            if N==16
                plot(r(:,1),r(:,2),'r-','LineWidth',3) %Plotning af løsning
            else
                plot(r(:,1),r(:,2),'b') %Plotning af løsning
            end
        end  
    end
    
r0=pos2(t0)+.5*b1*[0; -1; 0]'; %Startbetingelse
[s,r]=ode23(@difflign,sspan,r0,options); %L�sning af differentialligning
%plot(r0(1),r0(2),'kx')
plot(r(:,1),r(:,2)) %Plotning af l�sning
    
    
plot(0,-b1,'kx')
part2=pos2(d/(Beta*c));
plot(d,part2(2)+b2,'kx')
if -T0~=t0
    x=[];y=[];
    t=linspace(-T0,t0,100);
    for tt=1:100
        vec=pos2(t(tt));
        x=[x vec(1)]; y=[y vec(2)];
    end
    plot(x,y,'k--')
end
%axis equal
axis([-2e-8 2e-8 -2e-8 2e-8])
axis off
Film(ti)=getframe;
end


%% Prepare the new file.
    vidObj = VideoWriter('fieldlines.avi');
        vidObj.FrameRate = 4;

    open(vidObj);
 
    for k = 1:nframes
       % Write each frame to the file.
       writeVideo(vidObj,Film(k));
    end
    % Close the file.
    close(vidObj);
 
%%Differentialligning for feltlinierne
function drds=difflign(s,r)
% Beregner først den retarderede tid
global  e1 c gamma t0
tret=@(tr)sqrt(sum((r(1:3)'-pos2(tr)).^2))-c*(t0-tr);
options2=optimset('TolX',1E-24,'TolFun',1E-15);
TR=fzero(tret,0,options2);

% Beregner normalvektoren fra partikens retarderede position til
% feltpunktet
normalunorm=(r(1:3)'-pos2(TR));
norm=sqrt(normalunorm*normalunorm');
normal=normalunorm/norm;

%Beregner det elektriske felt
betavek=beta2(TR);
denom=(1-betavek*normal')^3;

Efelt=e1*((normal-betavek)/((gamma*(t0-TR)*c)^2*denom))  +  ...
    e1/c*(cross(normal, cross(normal-betavek,dbeta2(TR)))/(denom*c*(t0-TR)));
Enorm=sqrt(Efelt*Efelt');

drds=Efelt'/Enorm;



