close all; clear all;

global q B

q = 1;
B = [0; 0; 0.5];
m = 0.1;

v = 4*[10; 0; 0]; v2 = v;

x = -1; x2 = x;
y = -0.5; y2 = y + 0.01;

dt = 0.001;

r = m*v / (q*norm(B))

f = figure;
h1 = animatedline('color','blue','markersize',25,'marker','.');
h2 = animatedline('color','red','markersize',25,'marker','.');

h11 = animatedline('color','blue','linestyle','-','linewidth',2.5);
h22 = animatedline('color','red','linestyle','-','linewidth',2.5);

axis tight manual 
axis([-1, 1, -1 ,0]);
hold on

vid = VideoWriter('../../presentation/figures/deflection2');

open(vid);

plot([0.5, 0.5],[-1, 0],'k--','linewidth',1.5)
plot([-0.5, -0.5],[-1, 0],'k--','linewidth',1.5)
fill([-0.5, 0.5, 0.5, -0.5],[-1, -1, 0, 0],'b','FaceAlpha',0.2,'edgecolor','none')
text(-0.95,-0.1,'B = 0','FontSize',22)
text(0.55,-0.1,'B = 0','FontSize',22)
text(-0.25,-0.1,'B > 0','FontSize',22)

axis off

while x < 1 && y > -1
   clearpoints(h1);
   clearpoints(h2);
   
   if x < 0.5 && x > -0.5; a = 1/m * magneticforce(v); else; a = [0;0;0]; end
   if x2 < 0.5 && x2 > -0.5; a2 = -1/m * magneticforce(v2); else; a2 = [0;0;0]; end

   dv = a * dt;
   v = v + dv;
   dr = v * dt;
   x = x + dr(1);
   y = y + dr(2);

   dv2 = a2 * dt;
   v2 = v2 + dv2;
   dr2 = v2 * dt;
   x2 = x2 + dr2(1);
   y2 = y2 + dr2(2);

   addpoints(h1, x, y);
   addpoints(h11, x, y);
   addpoints(h2, x2, y2);
   addpoints(h22, x2, y2);
  
   drawnow
   
  frame = getframe(f);
  writeVideo(vid,frame);
  
end

close(vid);

function F = magneticforce(v)

global q B

F = q * cross(v, B);

end