clear; clc; close all;
format long

%% Calculating r
rE = 6378.137; mu=398600.4418;

a1 = 13000; e1 = 0.3;
a2 = 7226.58; e2 = 0.444819;
dw = deg2rad(301.901-50);

p1 = a1 * (1 - e1^2); h1 = sqrt(mu*p1);
p2 = a2 * (1 - e2^2); h2 = sqrt(mu*p2);

alpha=e2*cos(dw)-e1;
beta=e1*p2-e2*p1*cos(dw);
gamma=e1*e2*sin(dw);

A = 1 - 1/e1^2 - (alpha/gamma)^2;
B = 2*p1/e1^2 - 2*alpha*beta/gamma^2;
C = -((p1/e1)^2 + (beta/gamma)^2);

r = (-B+sqrt(B^2 - 4*A*C))/(2*A);
% r = (-B-sqrt(B^2 - 4*A*C))/(2*A);


SIN = (alpha*r+beta)/(gamma*r);
COS = (1/e1)*(p1/r - 1);
f = asin(SIN);
% acosd(COS);

aop1 = deg2rad(50);
aop2 = deg2rad(301.901);
I = deg2rad(20);
raan = deg2rad(30);

r_peri = r*[cos(f) sin(f) 0]';
r_eci = inv(R3(aop1)*R1(I)*R3(raan))*r_peri;


% r_peri = r*[cos(f-dw) sin(f-dw) 0]';
% r_eci = inv(R3(aop2)*R1(I)*R3(raan))*r_peri;

vr1 = (h1/r)*(e1*sin(f)/(1+e1*cos(f)));
vr2 = (h2/r)*(e2*sin(f-dw)/(1+e2*cos(f-dw)));
vt1 = h1/r;
vt2 = h2/r;

v_orb1 = [vr1 vt1 0]'; v1 = norm(v_orb1);
v_orb2 = [vr2 vt2 0]'; v2 = norm(v_orb2);

v_eci1 = inv(R3(aop1+f)*R1(I)*R3(raan))*v_orb1;
v_eci2 = inv(R3(aop1+f)*R1(I)*R3(raan))*v_orb2;

norm(v_eci1-v_eci2)










%% Graphing
% f_graph = 1:0.1:360;
% for i=1:size(f_graph,2)
%     radius1(i) = p1/(1+e1*cos(f_graph(i)));
%     radius2(i) = p2/(1+e2*cos(f_graph(i)-dw));
% end
% x1 = radius1.*cos(f_graph); y1 = radius1.*sin(f_graph);
% x2 = radius2.*cos(f_graph); y2 = radius2.*sin(f_graph);
% 
% 
% vx1 = cos(f)*vr1 - sin(f)*vt1;
% vy1 = sin(f)*vr1 + cos(f)*vt1;
% vx2 = cos(f)*vr2 - sin(f)*vt2;
% vy2 = sin(f)*vr2 + cos(f)*vt2;
% 
% figure;
% plot(x1,y1,'color','green'); hold on;
% plot(x2,y2,'color','#cc66ff'); hold on;
% % quiver(0,0,r_peri1(1),r_peri1(2),'color','green');hold on;
% % quiver(r_peri2(1),r_peri2(2),vr1,vt1)
% plot([0 r_peri(1)],[0 r_peri(2)],'color','red');
% plot([r_peri(1) r_peri(1)+500*vx1],[r_peri(2) r_peri(2)+500*vy1],'color','green');
% plot([r_peri(1) r_peri(1)+500*vx2],[r_peri(2) r_peri(2)+500*vy2],'color','#cc66ff');





function R1 = R1(angle)
    R1 = [ 1 0 0
        0 cos(angle) sin(angle)
        0 -sin(angle) cos(angle)];
end
function R3 = R3(angle)
    R3 = [ cos(angle) sin(angle) 0
        -sin(angle) cos(angle) 0
        0 0 1];
end


