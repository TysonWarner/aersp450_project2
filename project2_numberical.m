clear, clc, close all

load IODMeasurements2.mat
format longG

for i = 1 %1:3:22
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
    v2 = gibbs(rr1, rr2, rr3);
end


function v2 = gibbs(rr1, rr2, rr3)
    validate = dot(rr1/norm(rr1),cross(rr2,rr3)/norm(cross(rr2,rr3)));
    if ~(validate < 10^-8 && validate > -10^-8)
        disp('Error in Gibbs');
    end

    mu = 3.986*10^5;

    n = norm(rr1)*cross(rr2,rr3) + norm(rr2)*cross(rr3,rr2) + norm(rr3)*cross(rr1,rr2);
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
    for i = 1:length(r2_magnitude)
        j=1;
        if isreal(r2_magnitude(i)) && r2_magnitude(i) > 0
            r2(j) = r2_magnitude(i);
            j = j+1;
        end
    end
    %r2 = r2_magnitude(4); % positive and real

    rho1 = 1/D0*((6*(D31*tau1/tau3+D21*tau13/tau3)*r2^3 + mu*D31*(tau13^2-tau1^2)*tau1/tau3)/(6*r2^3+mu*(tau13^2-tau3^2))-D11);
    rho2 = A + mu*B*r2^(-3);
    rho3 = 1/D0*((6*(D13*tau3/tau1-D23*tau13/tau1)*r2^3 + mu*D13*(tau13^2-tau3^2)*tau3/tau1)/(6*r2^3+mu*(tau13^2-tau1^2))-D33);

    rr1 = R1 + double(rho1).*L1;
    rr2 = R2 + double(rho2).*L2;
    rr3 = R3 + double(rho3).*L3;
end


% function totalSeconds = getTime(timeStr)
%     timeComponents = textscan(timeStr, '%d:%d:%d:%f %s');
% 
%     days = timeComponents{1};
%     hours = timeComponents{2};
%     minutes = timeComponents{3};
%     seconds = timeComponents{4};
%     ampm = char(timeComponents{5});
% 
%     if strcmp(ampm, 'PM') && hours < 12
%         hours = hours + 12;
%     elseif strcmp(ampm, 'AM') && hours == 12
%         hours = 0;
%     end
% 
%     totalSeconds = days * 86400 + hours * 3600 + minutes * 60 + seconds;
% end
