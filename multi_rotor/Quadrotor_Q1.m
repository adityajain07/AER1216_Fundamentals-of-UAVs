%% Author: Aditya Jain
%  Date  : 8th December, 2021
%  About : AER1216 Project, Quadrotor

clc
clear
close all

m      = 0.42;            % mass of the quadrotor
g      = 9.8;             % gravity constant
W      = m*g;             % weight of the quadrotor
CD     = 0.97;   
S      = 0.01;            % reference area
pho    = 1.225;           % air density at sea level
D      = 8*0.0254;        % diameter of prop in meters
r_area = 4*(pi*(D/2)^2);  % total rotor area

%% Hovering based on APC propeller measurements
% The thrust should be equal to about 1.03 N for every propeller (4.12/4)
% CT and CP were obtained from the UIUC database for T = 1.07
% T   = CT*pho*n^2*D^4        % thrust for 1 prop
% RPM = 3581;
% n   = RPM/60;

CT  = 0.1449;
CP  = 0.0822;
T   = W/4;         % thrust required per propeller
n   = sqrt(T/(CT*pho*D^4));
RPM = n*60;
P   = CP*pho*n^3*D^5;        % power for 1 prop

P_quad = 4*P;                % power for 4 prop
X      = ['Power required to hover considering APC propellers measured characteristics: ', num2str(P_quad), ' W'];
disp(X)

Eb     = 3*3.7*1.5*3600;  % energy in the battery
eta_m  = 0.75;            % motor efficiency
eta_e  = 0.85;            % esc efficiency
te     = (Eb*eta_m*eta_e)/P_quad;
te_min = te/60;

X = ['Hover endurance from APC measurment (0th battery model) is: ', num2str(te_min), ' min'];
disp(X)

%% Endurance and Range based for Forward Speeds

V_for   = (0:0.5:40);          % forward speeds to check for
T       = zeros(size(V_for));  % thrust required
P_ind   = zeros(size(V_for));  % induced power
P_tot   = zeros(size(V_for));  % total power
P_tot_v = zeros(size(V_for));  % total power per velocity

for i = 1:length(V_for)
   V       = V_for(i);
   D       = 0.5*pho*S*CD*V^2;
   alpha_d = atan(D/W);
   T(i)    = sqrt(W^2 + D^2);
   
   % calculation for induced velocity
   p     = [1 2*V*sin(alpha_d) V^2 0 -(W^2 + D^2)/(2*pho*r_area)^2];
   rts   = roots(p);
   
   for j=1:length(rts)
      flag = isreal(rts(j));
      if flag && rts(j)>0
        v_ind = rts(j);
      end
   end
   P_ind(i)  = T(i)*v_ind;
   P_tot(i)  = T(i)*(v_ind + V*sin(alpha_d));
   P_tot_v(i) = P_tot(i)/V;
   
end

figure
plot(V_for, P_tot)
xlabel('Velocity (m/sec)')
ylabel('Power (W)')
title('Total Power v/s Velocity')
hold on 
plot(V_for(15), P_tot(15), 'ro', 'MarkerSize', 10)
text(V_for(15), P_tot(15),'Min. Total Power')
hold off

figure
plot(V_for, P_tot_v)
xlabel('Velocity (m/sec)')
ylabel('Power/Velocity (sec*W/m)')
title('Total Power/Velocity v/s Velocity')
hold on 
plot(V_for(20), P_tot_v(20), 'ro', 'MarkerSize', 10)
text(V_for(20), P_tot_v(20),'Min. Total Power/Velocity')
hold off

X = ['The required forward speed for max endurance is: ', num2str(V_for(15)), ' m/sec'];
disp(X)

X = ['The required forward speed for max range is: ', num2str(V_for(20)), ' m/sec'];
disp(X)

% Endurance Calculation
Eb        = 3*3.7*1.5*3600;  % energy in the battery
eta_m     = 0.75;            % motor efficiency
eta_e     = 0.85;            % esc efficiency
Pmin_end  = P_tot(15);       % min. power for endurance
te        = (Eb*eta_m*eta_e)/Pmin_end;
te_min    = te/60;

X = ['The maximum flying endurance is: ', num2str(te_min), ' min'];
disp(X)

% Range Calculation
Pmin_ran  = P_tot_v(20)*V_for(20);
te        = (Eb*eta_m*eta_e)/Pmin_ran;
ran_max   = te*V_for(20);
X         = ['The maximum range is: ', num2str(ran_max/1000), ' km'];
disp(X)

