clear; clc; close all;
format short

%% Calculating r
rE = 6378.137;

a1 = 13000; e1 = 0.3;
a2 = 7226.58; e2 = 0.444819;
dw = deg2rad(301.901-50);

p1 = a1 * (1 - e1^2);
p2 = a2 * (1 - e2^2);

alpha=e2*cos(dw)-e1;
beta=e1*p2-e2*p1*cos(dw);
gamma=e1*e2*sin(dw);

A = 1 - 1/e1^2 - (alpha/gamma)^2;
B = 2*p1/e1^2 - 2*alpha*beta/gamma^2;
C = -((p1/e1)^2 + (beta/gamma)^2);

r = (-B+sqrt(B^2 - 4*A*C))/(2*A);
% r = (-B-sqrt(B^2 - 4*A*C))/(2*A);


SIN = (alpha*r+beta)/(gamma*r)
COS = (p2/p1 + e2*SIN*sin(dw) - 1)/(e2*cos(dw) - (p2*e1)/p1)
asind(SIN)
acosd(COS)

%% Calculating position and velocity vectors
r_orbital = r*[1 0 0]';

syms theta I Omega omega f
A = R3(theta+f)*R1(I)*R3(Omega);
r_eci = simplify(inv(A)*r_orbital);
r_eci = subs(r_eci,I,deg2rad(20));
r_eci = subs(r_eci,Omega,deg2rad(30));












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


