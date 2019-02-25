clear all; close all;

data = load("../src/test.txt");
x = data(:,1);
y = data(:,2);
z = data(:,3);

figure
hold on
plot(x(z == 5), y(z == 5),'x')
axis equal
box on
grid on
