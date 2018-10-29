close all; clear all; clc

Npoints = 100;
lb = -4; ub = 4;
assert(lb < ub)

r = linspace(lb, ub, Npoints);
[X,Y] = meshgrid(r);
centroid = 0;
d = 3.57;

center = [0, d; d, d; -d, d; d, -d; -d, -d; d, 0; 0, 0; -d, 0; 0, -d]

UU = zeros(size(X));
for i = 1:length(center)
    
   xcenter = center(i, 1);
   ycenter = center(i, 2);
   field = strcat('R',num2str(i));
   R.(field) = sqrt((X - xcenter).^2 + (Y - ycenter).^2);
   U.(field) = potential(R.(field));
   UU = UU + U.(field).DT;
   
end

r = linspace(lb, ub, 100*Npoints);
u1 = potential(r,-2);
u2 = potential(r,2);
E0 = e_field(r,0);
E1 = e_field(r,d);
E2 = e_field(r,-d);

f = figure;
subplot(1,2,1)
[C, h] = contour(X, Y, -UU, 1000);
grid on
xlabel('x pos [\AA]','fontsize',22,'interpreter','latex')
ylabel('y pos [\AA]','fontsize',22,'interpreter','latex')
c = colorbar();
colormap(flipud(parula))
c.Label.String = 'U [eV]';
c.Label.FontSize = 22;
c.Label.Interpreter = 'latex';
set(gca, 'FontSize', 18)
box on
title('a) axial potential','fontsize',22,'interpreter','latex')
% set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
% f = figure;
subplot(1,2,2)
title('b) axial field','fontsize',22,'interpreter','latex')
hold on
plot(r, E0 + E1 + E2,'-','linewidth',1.5)
xlabel('Transverse distance [\AA]','fontsize',22,'interpreter','latex')
ylabel('E [V/\AA]','fontsize',22,'interpreter','latex')
ax = gca;
ax.YAxisLocation = 'right';
box on
grid on
set(gca, 'FontSize', 18)
set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[36, 12],'PaperPosition',[0, 0, 36, 12],'Position',[0 0 36 12])
print(f,'../../figures/axial_potential_field.pdf', '-dpdf','-r600','-painters')

function U = potential(r, centroid)
    Z1 = 1; % electron/positron
    Z2 = 6; % carbon
    d = 3.56; % interatomic spacing <100> C [Å]
    C = 3;
    a0 = 0.529177211; % bohr radius [Å]
    e2 = 14.4; % e^2 in  [eV Å] (cgs units)
    rho2 = 0.043^2; % mean squared thermal displacement of string [Å]
    aH = 0.89 * a0 * (Z1^(2/3) + Z2^(2/3))^(-1/2); % [Å]
    
    b = [36.9951,11.2966,2.8139,0.3456];
    a = [0.7307, 1.1951, 0.4563, 0.1247];

    B = b / (4*pi^2);

    if min(size(r)) == 1
        r = r - centroid;
    end

    U1 = 0;
    for i = 1:length(B)

        U1 = U1 + a(i)/(B(i) + rho2) * exp(- r.^2 / (B(i) + rho2));

    end

    U.DT = 2*Z1*e2*a0/d * U1;
    U.H = (Z1*Z2*e2)/d * log(1 + C^2 * aH^2 ./ (r.^2 + 0.5*rho2));
end

function E = e_field(r, centroid)
    Z1 = 1; % electron/positron
    a0 = 0.529177211; % bohr radius [Å]
    e2 = 14.4; % e^2 in  [eV Å] (cgs units)
    rho2 = 0.043^2; % mean squared thermal displacement of string [Å]
    d = 3.56; % interatomic spacing <100> C [Å]
    b = [36.9951,11.2966,2.8139,0.3456];
    a = [0.7307, 1.1951, 0.4563, 0.1247];

    B = b / (4*pi^2);

    if min(size(r)) == 1
        r = r - centroid;
    end
    
    E1 = 0;
    for i = 1:length(B)

        E1 = E1 + a(i)/(B(i) + rho2)^2 * r .* exp(- r.^2 / (B(i) + rho2));

    end
    
    E = -(-4*Z1*e2*a0)/d * E1;
    
end
