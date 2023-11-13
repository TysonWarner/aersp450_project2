clear, clc, close all

load IODMeasurements2.mat
format longG

mu = 3.986004418*10^5;

plotEarth; hold on;

for i = 1:3:22 %19:3:22
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

    if (i == 7 || i == 13 || i == 22)
        delete(graph)
    end
end

legend('Earth''s Surface','Orbit 1', 'Orbit 2', 'Orbit 3', 'Orbit 4', 'Orbit 5')


function plotEarth
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
    mu = 3.986004418*10^5;
    state0 = [r0 v0];
    T = 2*pi*sqrt(a^3/mu);
    tspan = [0, T]; 
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
    [t, state] = ode78(@TwoBP, tspan, state0, options);

    % graph = plot3(state(:,1),state(:,2),state(:,3));
    graph = quiver3(state(:,1),state(:,2),state(:,3),state(:,4),state(:,5),state(:,6),5);
    % plot3(state(:,1), state(:,2), state(:,3), 'color', 'red');
    title('Orbit in ECI Frame');
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
