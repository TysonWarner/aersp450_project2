%% 
%  1) Initial calculations found at beginning of the main loop.
%     Gauss algorithm located in the Gauss function.
%  2) The results from the Gauss alogrithm are found in the main report as
%     well as the reasoning.
%  3) Gibb's method is found in the function Gibb. The resulting velocities
%     can be found in the report.
%  4) The conversion from r2 and v2 to orbital elements is found in the 
%     RV2OE function. These values can be found in the main report.
%  6) The orbits are plotting using the graph function. The two-body
%     problem was solved to accomplish this. The resulting plot can be found
%     in the main report.
%  7) The orbital elements were used to calculate the delta-Vs. The code
%     was derived in the main report and can be found implemented in the
%     deltaVcalc function. The results can be found in the main report.


%% main
clear, clc, close all

load IODMeasurements2.mat
format longG

mu = 3.986*10^5;
datasetCt = 1;

surfaceEarth; hold on;

for i = 1:3:22 
    L1 = [cosd(AZIMUTH(i))*cosd(ELEVATION(i)) sind(AZIMUTH(i))*cosd(ELEVATION(i)) sind(ELEVATION(i))];
    L2 = [cosd(AZIMUTH(i+1))*cosd(ELEVATION(i+1)) sind(AZIMUTH(i+1))*cosd(ELEVATION(i+1)) sind(ELEVATION(i+1))];
    L3 = [cosd(AZIMUTH(i+2))*cosd(ELEVATION(i+2)) sind(AZIMUTH(i+2))*cosd(ELEVATION(i+2)) sind(ELEVATION(i+2))];
    R1 = RSITES(i,:);
    R2 = RSITES(i+1,:);
    R3 = RSITES(i+2,:);
    t1 = TIMES(i);
    t2 = TIMES(i+1);
    t3 = TIMES(i+2);

    [rr1, rr2, rr3] = gauss(t1, t2, t3, L1, L2, L3, R1, R2, R3);
    
    if size(rr1,1) == 3 %if 3 possible orbits from gauss
        if i == 19 % One set of data where this occurs and the correct v2 is selected
            j = 3;
            v2(j,1:3) = gibbs(rr1(j,1:3), rr2(j,1:3), rr3(j,1:3));
            [a, eNORM, I, RAAN, AOP, f] = RV2OE(rr2(j,1:3),v2(j,1:3),mu);
            graph = plotOrbit(rr2(j,1:3),v2(j,1:3),a);
        elseif i == 22 % One set of data where this occurs and the correct v2 is selected
            j = 2;
            v2(j,1:3) = gibbs(rr1(j,1:3), rr2(j,1:3), rr3(j,1:3));
            [a, eNORM, I, RAAN, AOP, f] = RV2OE(rr2(j,1:3),v2(j,1:3),mu);
            graph = plotOrbit(rr2(j,1:3),v2(j,1:3),a);
        % else 
        %     for j = 2 % plots were compared for the cases aboved and selected.
        %         v2(j,1:3) = gibbs(rr1(j,1:3), rr2(j,1:3), rr3(j,1:3));
        %         [a, eNORM, I, RAAN, AOP, f] = RV2OE(rr2(j,1:3),v2(j,1:3),mu);
        %         plotOrbit(rr2(j,1:3),v2(j,1:3),a);
        %     end
        end
    elseif size(rr1,1) == 1
        v2 = gibbs(rr1, rr2, rr3);
        [a, eNORM, I, RAAN, AOP, f] = RV2OE(rr2,v2,mu);
        graph = plotOrbit(rr2,v2,a);
    end

    orbElemsHist(datasetCt,1:6) = [a, eNORM, I, RAAN, AOP, f];

    datasetCt = datasetCt+1;

    if (i == 7 || i == 13 || i == 22)
        delete(graph)
    end
end

legend('Earth''s Surface','Orbit 1', 'Orbit 2', 'Orbit 3', 'Orbit 4', 'Orbit 5');

disp('Refer to excel document or main report document for tables of values.');

%% Calculating the delta-Vs
for n = 2:8
    a1 = orbElemsHist(n-1,1);
    e1 = orbElemsHist(n-1,2);
    I = orbElemsHist(n-1,3);
    raan = orbElemsHist(n-1,4);
    aop1 = orbElemsHist(n-1,5);

    a2 = orbElemsHist(n,1);
    e2 = orbElemsHist(n,2);
    aop2 = orbElemsHist(n,5);

    deltaVlist(n-1,1) = deltaVcalc(a1, e1, aop1, a2, e2, aop2, I, raan);
end

%% Main Functions
function deltaV = deltaVcalc(a1, e1, aop1, a2, e2, aop2, I, raan)
    mu = 3.986*10^5;

    dw = aop2-aop1;

    p1 = a1 * (1 - e1^2); h1 = sqrt(mu*p1);
    p2 = a2 * (1 - e2^2); h2 = sqrt(mu*p2);

    alpha=e2*cos(dw)-e1;
    beta=e1*p2-e2*p1*cos(dw);
    gamma=e1*e2*sin(dw);

    A = 1 - 1/e1^2 - (alpha/gamma)^2;
    B = 2*p1/e1^2 - 2*alpha*beta/gamma^2;
    C = -((p1/e1)^2 + (beta/gamma)^2);

    r = (-B+sqrt(B^2 - 4*A*C))/(2*A);

    SIN = (alpha*r+beta)/(gamma*r);
    f = asin(SIN);

    vr1 = (h1/r)*(e1*sin(f)/(1+e1*cos(f)));
    vr2 = (h2/r)*(e2*sin(f-dw)/(1+e2*cos(f-dw)));
    vt1 = h1/r;
    vt2 = h2/r;

    v_orb1 = [vr1 vt1 0]'; 
    v_orb2 = [vr2 vt2 0]'; 

    v_eci1 = inv(R3ang(aop1+f)*R1ang(I)*R3ang(raan))*v_orb1;
    v_eci2 = inv(R3ang(aop1+f)*R1ang(I)*R3ang(raan))*v_orb2;

    deltaV = norm(v_eci1-v_eci2);
end

function surfaceEarth
    % method to plot Earth
    Re = 6378;
    [Xe, Ye, Ze] = sphere;
    Xe = Xe * Re;
    Ye = Ye * Re;
    Ze = Ze * Re;

    %figure;
    surf(Xe, Ye, Ze, 'FaceColor', 'b');
end

function graph = plotOrbit(r0,v0,a)
    % Propogate orbit using ode45 in ECI frame
    mu = 3.986*10^5;
    state0 = [r0 v0];
    T = 2*pi*sqrt(a^3/mu);
    tspan = [0, T]; 
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t, state] = ode78(@TwoBP, tspan, state0, options);

    %plot3(state(:,1),state(:,2),state(:,3));
    graph = plot3(state(:,1), state(:,2), state(:,3));
    %quiver3(state(:,1),state(:,2),state(:,3),state(:,4),state(:,5),state(:,6));

    %title('Orbit in ECI Frame');
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    grid on;
    % xlim([-3e4, 3e4]);
    % ylim([-3e4, 3e4]);
    % zlim([-3e4, 3e4]);
    axis equal;
end

function statedot = TwoBP(t, state)
    mu = 3.986*10^5;
    
    r = state(1:3);
    v = state(4:6);

    rdot = v;
    rdoubledot = -mu * r / norm(r)^3;
    
    statedot = [rdot; rdoubledot];
end

function v2 = gibbs(rr1, rr2, rr3)
    validate = dot(rr1/norm(rr1),cross(rr2,rr3)/norm(cross(rr2,rr3)));
    if ~(validate < 10^-8 && validate > -10^-8)
        disp('Error in Gibbs');
    end

    mu = 3.986*10^5;

    n = norm(rr1)*cross(rr2,rr3) + norm(rr2)*cross(rr3,rr1) + norm(rr3)*cross(rr1,rr2);
    d = cross(rr1,rr2) + cross(rr2,rr3) + cross(rr3,rr1);
    s = rr1*(norm(rr2)-norm(rr3)) + rr2*(norm(rr3)-norm(rr1)) + rr3*(norm(rr1)-norm(rr2));

    v2 = sqrt(mu/dot(n,d))*(cross(d,rr2)/norm(rr2) + s);
end

function [rr1, rr2, rr3] = gauss(t1, t2, t3, L1, L2, L3, R1, R2, R3)
    mu = 3.986*10^5;

    tau1 = t1-t2;
    tau3 = t3-t2;
    tau13 = t3-t1;

    D0 = dot(L1,cross(L2,L3));
    D11 = dot(R1,cross(L2,L3));
    D12 = dot(R1,cross(L1,L3));
    D13 = dot(R1,cross(L1,L2));
    D21 = dot(R2,cross(L2,L3));
    D22 = dot(R2,cross(L1,L3));
    D23 = dot(R2,cross(L1,L2));
    D31 = dot(R3,cross(L2,L3));
    D32 = dot(R3,cross(L1,L3));
    D33 = dot(R3,cross(L1,L2));

    A = 1/D0*(-tau3*D12/tau13 + D22 + tau1/tau13*D32);
    B = 1/6/D0*(-(tau13^2-tau3^2)*tau3/tau13*D12 + (tau13^2-tau1^2)*tau1/tau13*D32);

    a = -A^2 - 2*A*dot(L2,R2)-norm(R2)^2;
    b = -2*mu*B*(A+dot(L2,R2));
    c = -mu^2*B^2;

    r2_magnitude = roots([1.0, 0, double(a), 0, 0, double(b), 0, 0, double(c)]);
    j=1;
    for i = 1:length(r2_magnitude)
        if isreal(r2_magnitude(i)) && r2_magnitude(i) > 0
            r2(j) = r2_magnitude(i);
            j = j+1;
        end
    end
    
    for k = 1:j-1
        rho1(k) = 1/D0*((6*(D31*tau1/tau3+D21*tau13/tau3)*r2(k)^3 + mu*D31*(tau13^2-tau1^2)*tau1/tau3)/(6*r2(k)^3+mu*(tau13^2-tau3^2))-D11);
        rho2(k) = A + mu*B*r2(k)^(-3);
        rho3(k) = 1/D0*((6*(D13*tau3/tau1-D23*tau13/tau1)*r2(k)^3 + mu*D13*(tau13^2-tau3^2)*tau3/tau1)/(6*r2(k)^3+mu*(tau13^2-tau1^2))-D33);

        rr1(k,1:3) = R1 + double(rho1(k)).*L1;
        rr2(k,1:3) = R2 + double(rho2(k)).*L2;
        rr3(k,1:3) = R3 + double(rho3(k)).*L3;
    end
end

function [a, eNORM, I, RAAN, AOP, f] = RV2OE(r,v,mu)
    % RV2OE - A function that accepts inputs position, velocity 
    % in ECI frame, and μ value, and outputs the orbital elements

    % Solving for semi-major axis
    energy = (norm(v)^2)/2 - mu/(norm(r));
    a = -mu/(2*energy);
    
    % Solving for eccentricity
    h = cross(r,v);
    p = (norm(h))^2/mu;
    e = (cross(v,h)/mu)-(r/norm(r));
    eNORM = norm(e);

    % Solving for true anomaly
    fpa = asin((dot(r,v))/(norm(r)*norm(v)));
    if fpa > 0
        f = acos((1/norm(e))*((p/norm(r))-1));
    else
        f = -acos((1/norm(e))*((p/norm(r))-1));
    end

    % Solving for inclination
    I = acos(h(3)/norm(h));

    % Solving for right ascension of the ascending node
    n = cross([0 0 1],h)/norm(cross([0 0 1],h));    % Node Vector
    sinRAAN = dot(n,[0 1 0]);
    cosRAAN = dot(n,[1 0 0]);
    if cosRAAN<0                                   % Finding RAAN
        RAAN = atan(sinRAAN/cosRAAN)+pi;
    else
        RAAN = atan(sinRAAN/cosRAAN);
    end

    % Solving for argument of perigee
    if dot(e,[0 0 1])<0  % Finding Argument of Periapsis
        AOP = -acos(dot(e,n)/norm(e));
    else
        AOP = acos(dot(e,n)/norm(e));
    end
end

%% Auxillary Functions
function R1ang = R1ang(angle)
    R1ang = [ 1 0 0
        0 cos(angle) sin(angle)
        0 -sin(angle) cos(angle)];
end

function R3ang = R3ang(angle)
    R3ang = [ cos(angle) sin(angle) 0
        -sin(angle) cos(angle) 0
        0 0 1];
end
