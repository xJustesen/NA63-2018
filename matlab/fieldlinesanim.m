clear all; clc; close all
%Definitioner
global e c m Beta gamma t0  dt F
e=1.6*2.99*1E-10; %elektron
c=2.99E10; % Lysets hastighed i cm/s
m=9.11E-28; %Elektronens masse
Beta = 0.95; %Beta faktoren
gamma=sqrt(1-Beta^2)^(-1); %Gamma faktoren
dt = 1e-16;
F = gamma * m * (Beta/dt + gamma^2 * Beta^3 / dt / c^2)

%% Løs af differentialligning for forskellige startbetingelser

sspan=[-1E-10 1E-10]; %Løsningsområde for difflign
options=odeset('Abstol',1E-15,'Reltol',3E-14); %Indstillinger for ode45
nframes = 60;
n=8; % Det halve antal feltlinier

f = figure;
hold on
set(gca,'PlotBoxAspectRatio',[1 1 1])
% axis([-2e-9 2e-9 -2.5e-10 2.5e-10])
axis off
filename = 'testAnimated';

for ti=1:nframes
    
    clf
    hold on
    set(gca,'PlotBoxAspectRatio',[1 1 1])
%     axis([-2e-9 2e-9 -2.5e-10 2.5e-10])

    t0 = 1e-18 * (1.5 * ti - 0.5 * nframes)/nframes; %Tidspunktet for billedet af feltlinierne
    
    for N = 0:n-1
        
        r0 = bane(t0, gamma) + 0.5 * [cos(atan(gamma * tan(pi * N / n))); sin(atan(gamma*tan(pi*N/n))); 0] * 0.2E-10; %Startbetingelse
        difflign2 = @(r, s) difflign(r, s, gamma);
        [~,r] = ode23(difflign2, sspan, r0, options); %Løsning af differentialligning
        plot(r(:,1),r(:,2),'k') %Plotning af løsning

        if N ~= n/2
            
        r0 = bane(t0, gamma) + 0.5 * [-cos(atan(gamma * tan(pi * N / n))); sin(atan(gamma*tan(pi*N/n))); 0] * 0.2E-10; %Startbetingelse
        difflign2 = @(r, s) difflign(r, s, gamma);
        [~,r] = ode23(difflign2, sspan, r0, options); %Løsning af differentialligning
            
            if N== floor((n-1)/2)
                plot(r(:,1),r(:,2),'r-','LineWidth',1.5) %Plotning af løsning
            else
                plot(r(:,1),r(:,2),'k') %Plotning af løsning
            end
            
        end
        
    end
    
    r0 = bane(t0, gamma) + 0.5 * [0; -1; 0] * 0.2E-10; %Startbetingelse
    difflign2 = @(r, s) difflign(r, s, gamma);
    [~,r] = ode23(difflign2, sspan, r0, options);
    
    plot(r(:,1),r(:,2),'k')
    
    % Gem plotvindue som billede
    frame = getframe(f);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Skriv til GIF fil
    if ti == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
    disp(['Frame no. ', num2str(ti),' of ', num2str(nframes), ' finished']);
%     drawnow

end


%% funktioner
function drds = difflign(s, r, gamma)
    global F e t0 c m;
    
    % Beregn retarderet tid
    tret=@(tr) sqrt(sum((r - bane(tr, gamma)).^2))-c*(t0-tr);
    options=optimset('TolX',1E-24,'TolFun',1E-15);
    TR = fzero(tret,t0,options);
    
    % Beregn position til retarderet tid
    [x, beta] = bane(TR, gamma);
    gamma = sqrt(1-norm(beta)^2)^(-1);
    dbeta = 1/gamma/m * ([F; 0; 0] - ([F, 0, 0] * beta) * beta);
    
    % Beregn normalvektoren fra partikens retarderede position til feltpunktet
    R = r - x;
    nhat = R / norm(R);
    
    % Beregn feltet
    denom = (1 - nhat' * beta)^3;
    Efelt = e * ((1 - norm(beta)^2)*(nhat - beta)/(norm(R)^2 * denom)  +  1/c * (cross(nhat, cross(nhat - beta, dbeta))) / (denom * norm(R)));

    % Returner diff. lign. for felt
    drds = Efelt / norm(Efelt);
end

function [x, beta] = bane(tret, gamma)
    global m c F Beta;
    
    % Beregn beta = v/c
    if tret > 0; beta = [Beta - F/m/gamma * tret/c; 0; 0]; else; beta = [Beta; 0; 0]; end
    if beta(1) < 0; beta(1) = 0; end
    
    % Beregn position x
    if tret > 0
        x = (-0.5 * 1/gamma/m * ([F; 0; 0] - ([F, 0, 0] * beta) * beta) * tret^2) + Beta*tret*c * [1; 0; 0];
    else
        x = Beta*tret*c * [1; 0; 0];
    end
    
end
