clear all; close all; clc;

beta = [0.1; 0];
t = linspace(0,1000,2000);
xmax = max(t) * beta(1);
posmax = max(t) * beta;
x = 0;
y = 0;

pos = [0, 0];

% h = animatedline('color','magenta','markersize',25,'marker','.');
% 
% axis([0,xmax,-1, 1]); % burde ikke hardcodes...
% box on

a = tic;
framerate = 60;

xsamples = linspace(-1, 1, 25);
ysamples = linspace(-1, 1, 25);
[X, Y] = meshgrid(xsamples, ysamples);

Ex = (1 ./ (X - pos(1)).^2);
Ey = (1 ./ (Y - pos(1)).^2);

Ex(Ex == +inf) = 1e3;
Ey(Ey == +inf) = 1e3;

quiver(X, Y, Ey, Ex,10)

% while any(pos < posmax)
%     % husk at ryde punkter hver opdatering så der ikke plottes oveni
%     clearpoints(h)
% 
%     if pos(1) > 0.5 * xmax
%        beta(2) = 0.001; 
%     end
%     
%     % opdater koordinat af punktladning   
%     pos = beta * mean(diff(t)) + pos;
% 
%     % opdater punkt til animation
%     addpoints(h, pos(1), pos(2));
% 
%     
%     % opdater animation
%     if toc(a) > 1/framerate
%         drawnow 
%         a = tic; % reset timer så vi opdaterer med konstant frekvens
%     end
% 
% end
% 
% drawnow; % drawnow efter plottet så scriptet automatisk stopper