close all; clear all; clc;

a = 1/10;
x = linspace(0, 50, 1000);

f = @(x) expo_inv_dist(a, x);
g = @(x) expo_dist(a, x);

u = rand(1,500);
X = inv_sample_analytical(f, u);
out = inv_sample_numeric(g, 500);

f = figure;
hold on
plot(X, u,'o')
plot(out(:,1), out(:,2),'x')
plot(x, expo_dist(a, x),'k--');
xlim([0, 50])
box on
grid on
legend({'Sampling (analytical)','Sampling (numeric)','f(x)'},'fontsize',18,'interpreter','latex','location','best')
ylabel('u','fontsize',22,'interpreter','latex');xlabel('x','fontsize',22,'interpreter','latex')
set(gca, 'fontsize',18)
set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
print(f,'../../figures/inv_trans_samp.pdf', '-dpdf','-r600','-painters')


function y = expo_inv_dist(a, x)
    y = -1/a * log(1-x);
end

function y = expo_dist(a, x)
    y = 1 - exp(-a*x);
end

function y = inv_sample_analytical(f, u)
    y = f(u);
end

function out = inv_sample_numeric(f, n)
    y = [];
    U = [];
    g = @(u, x) f(x) - u;
    u = rand(1, n);
    options = optimoptions('fsolve','Display','none');
    for i = 1:n
        x = fsolve(@(x) g(u(i), x), 0, options);
        y = [y; x];
        U = [U; u];
    end
    out = [y, u'];
end