close all; clear all; clc;

a = 1/10;
x = linspace(0, 50, 1000);

f = @(x) expo_inv_dist(a, x);
g = @(x) expo_dist(a, x);

u = rand(1,500);
X = inv_sample_analytical(f, u);
out = inv_sample_numeric(g, 500);
out2 = inv_sample_numeric(@gaussian_cdf,500);

f = figure;
hold on
plot(X, u,'o')
plot(out(:,1), out(:,2),'x')
plot(linspace(min(out(:,1)),max(out(:,1)),100), expo_dist(a, linspace(min(out(:,1)),max(out(:,1)),100)),'k--','linewidth',1.5);
xlim([0, 50])
box on
grid on
legend({'Sampling (analytical)','Sampling (numeric)','CDF(x)'},'fontsize',18,'interpreter','latex','location','best')
ylabel('u','fontsize',22,'interpreter','latex');xlabel('x','fontsize',22,'interpreter','latex')
set(gca, 'fontsize',18)
set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
print(f,'../../figures/inv_trans_samp.pdf', '-dpdf','-r600','-painters')

% f = figure;
% hold on
% % plot(X, u,'o')
% plot(out2(:,1), out2(:,2),'x')
% plot(linspace(min(out2(:,1)),max(out2(:,1)),100), gaussian_cdf(linspace(min(out2(:,1)),max(out2(:,1)),100)),'k--');
% % xlim([0, 50])
% box on
% grid on
% % legend({'Sampling (analytical)','Sampling (numeric)','f(x)'},'fontsize',18,'interpreter','latex','location','best')
% ylabel('u','fontsize',22,'interpreter','latex');xlabel('x','fontsize',22,'interpreter','latex')
% set(gca, 'fontsize',18)
% set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])


function y = gaussian_cdf(x)
    sigma = 1.0; mu = 0;
    y = 0.5 * (1 + erf( (x - mu)/(sigma * sqrt(2))));
end

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
    y = zeros(n,1);
    g = @(u, x) f(x) - u;
    u = rand(1, n);
    options = optimoptions('fsolve','Display','none');
    for i = 1:n
        x = fsolve(@(x) g(u(i), x), 0.5, options);
        y(i) = x;
    end
    out = [y, u'];
end