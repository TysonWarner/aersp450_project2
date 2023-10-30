clear; clc; close all;

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

r_same_1 = (-B+sqrt(B^2 - 4*A*C))/(2*A);
r_same_2 = (-B-sqrt(B^2 - 4*A*C))/(2*A);
r_a = a2 * (1 + e2);




f = deg2rad(0:.1:360);
r1 = zeros(size(f)); r2 = r1;

for i=1:size(f,2)
    r1(i) = p1 / (1 + e1 * cos(f(i)));
    r2(i) = p2 / (1 + e2 * cos(f(i)-dw));
end

plot(f,r1); hold on;
plot(f,r2); hold on;
yline([r_same_1 r_same_2])














