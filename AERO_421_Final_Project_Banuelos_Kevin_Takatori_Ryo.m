%% Final Project
% AERO 421
% Kevin Banuelos Ryo Takatori

% t_ means tumble
% dt_ means de-tumble
% no_ means normal operation

%% Assumptions
% Thruster fuel mass change is negligent
% Mass of thruster selected is within total mass
% Thruster doen't alter spacecraft geometry
% Magnetic dipole is same as HW 6&7
% Coefficient of drag is constant 2.2
% Air density based on NRLMSISE00
function AERO_421_Final_Project_Banuelos_Kevin_Takatori_Ryo
clc, clear all , close all
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Ode settings
%% Given
% Classical orbital elements
zp = 600; %
re = 6378;
mu = 398600; % Gravitational constant [km^3/s^2]
h = sqrt(mu*(zp+re)); % [km^2/s]
e = 0; % Eccentricity
OMEGA = 0; % RAAN [deg]
inc = 40; % Inclination [deg]
omega = 0; % Argument of perigee [deg]
theta = 0; % True anomaly [deg]
period = 2*pi*sqrt((zp+re)^3/mu);
t_w_bG_i = [0.05; -0.02; -0.03]; % Intial angular velocity [rad/s]
t_eps_bLVLH_i = [0;0;0]; % Intial quaternion vector
t_eta_bLVLH_i = 1; % Intial quaternion scalar
w_sc = 2; % Width of spacecraft [m]
h_sc = 2; % Height of spacecraft [m]
d_sc = 2; % Depth of spacecraft [m]
w_sensor = 0.25; % Width of sensor [m]
h_sensor = 0.25; % Height of sensor [m]
d_sensor = 1; % Depth of sensor [m]
m_sc_tot = 100; % Spacecraft mass [kg]

% Post Deployment/normal operation
no_w_bG_sc_i = [0.001; -0.001; 0.002];
m_sensor = m_sc_tot*0.1; % Sensor mass [kg]
m_solar = m_sc_tot*0.1; % Solar panel mass [kg]
m_cube = m_sc_tot*0.7-4.85*3;
w_solar = 3; % Width of solar panel [m]
h_solar = 2; % Height of solar panel [m]
d_solar = 0.05; % Depth of solar panel [m]
a_x = 2*2 + 0.25*1; % X-face area [m^2]
a_y = 2*2 + 0.25*1; % Y-face area [m^2]
a_z = 2*2 + 2*2*3; % Z-face area [m^2]
a_solar = 2*3;
a_cube = 2*2;
a_sensor_side = 1*0.25;
a_sensor_bottom = 0.25^2;
mdp = [0; 0; -0.5]; % Residual magnetic dipole [A*m^2]
C_D = 2.2; % Coefficient of drag
sp = 4.5*10^-6; % Solar Pressure @ 1 [AU] [N/m^2]
sun_ECI = [1;0;0]; % Sun vector

%% Tumbling
% Let the spacecraft tumble for one orbit to dealign spacecraft
% This will give the intial quaternion for de-tumble
% Moment of inertia of spacecraft in body frame
I_sc = (1/12)*m_sc_tot*[h_sc^2+d_sc^2 0 0;
    0 w_sc^2+d_sc^2 0;
    0 0 w_sc^2+h_sc^2]; % Moment of inertia for spacecraft
[ rperi, vperi, t_reci, t_veci ] = coe2rv( mu, e, h, inc, OMEGA, omega, theta); % Getting position and velocity from coes
% LVLH to Body transformation
t_C_bLVLH = quat2dcm([t_eta_bLVLH_i;t_eps_bLVLH_i]'); % Quaternion to DCM

% ECI to LVLH transformation
t_r = norm(t_reci); % Position in ECI
t_z_hat = -t_reci/t_r; % Z LVLH
t_y_hat = -cross(t_reci,t_veci)/norm(cross(t_reci,t_veci)); % Y LVLH
t_x_hat = cross(t_y_hat,t_z_hat); % X LVLH
t_C_LVLHG = [t_x_hat,t_y_hat,t_z_hat]'; % ECI to LVLH DCM

% ECI to body transformation
t_C_bG = t_C_bLVLH*t_C_LVLHG; % ECI to body DCM
t_eta_bG = sqrt(trace(t_C_bG)+1)/2; % Quaternion scalar

% Quaternion vectors
t_eps_bG_1 = (t_C_bG(2,3)-t_C_bG(3,2))/(4*t_eta_bG);
t_eps_bG_2 = (t_C_bG(3,1)-t_C_bG(1,3))/(4*t_eta_bG);
t_eps_bG_3 = (t_C_bG(1,2)-t_C_bG(2,1))/(4*t_eta_bG);
t_eps_bG = [t_eps_bG_1; t_eps_bG_2; t_eps_bG_3];

tinitial = 0; % Initial time [s]
tfinal = period; % Orbit time [s]
t_t = [tinitial tfinal]; % Time span [s]
state = [t_reci;t_veci; t_eps_bG; t_eta_bG; t_eps_bLVLH_i; t_eta_bLVLH_i;t_w_bG_i]; % State vectors
[t_t,t_statenew]=ode45(@tumble,t_t,state,options,mu,I_sc); % Ode45
t_reci_new = [t_statenew(:,1),t_statenew(:,2),t_statenew(:,3)]; % Position states [km]
t_veci_new = [t_statenew(:,4),t_statenew(:,5),t_statenew(:,6)]; % Velocity states [km/s]
% ECI to body
t_eps_bG_new = [t_statenew(:,7),t_statenew(:,8),t_statenew(:,9)]; % Quaternion vector states
t_eta_bG_new = t_statenew(:,10); % Quaternion scalar states
% LVLH to body
t_eps_bLVLH_new = [t_statenew(:,11),t_statenew(:,12),t_statenew(:,13)]; % Quaternion vector states
t_eta_bLVLH_new = t_statenew(:,14); % Quaternion scalar states
t_w_new = [t_statenew(:,15),t_statenew(:,16),t_statenew(:,17)]; % Angular velocity states [rad/s]

% Euler angles for ECI to body
k = 0;
j = 0;
for i = 1:length(t_t)
    C_bG_new = quat2dcm([t_eta_bG_new(i),t_eps_bG_new(i,:)]);
    [t_yaw_bG_new(i),t_pitch_bG_new(i),t_roll_bG_new(i)] = dcm_to_ypr(C_bG_new);

    %Angle fixing beyond 180 due matlab atan2d restrictions
    if i>1
        if t_yaw_bG_new(i)-(t_yaw_bG_new(i-1)+k*360)>170
            k = k+1;
        end
        if (t_yaw_bG_new(i-1)+k*360)-t_yaw_bG_new(i)>170
            k = k-1;
        end
        t_yaw_bG_new(i) = t_yaw_bG_new(i) - k*360;        
        if t_roll_bG_new(i)-(t_roll_bG_new(i-1)+j*360)>170
            j = j+1;
        end
        if (t_roll_bG_new(i-1)+j*360)-t_roll_bG_new(i)>170
            j = j-1;
        end
        t_roll_bG_new(i) = t_roll_bG_new(i) -j*360;
    end
end
% Euler angles for LVLH to body
k = 0;
j = 0;
for i = 1:length(t_t)
    C_bLVLH_new = quat2dcm([t_eta_bLVLH_new(i),t_eps_bLVLH_new(i,:)]);
    [t_yaw_bLVLH_new(i),t_pitch_bLVLH_new(i),t_roll_bLVLH_new(i)] = dcm_to_ypr(C_bLVLH_new);

    %Angle fixing beyond 180 due matlab atan2d restrictions
    if i>1
        if t_yaw_bLVLH_new(i)-(t_yaw_bLVLH_new(i-1)+k*360)>170
            k = k+1;
        end
        if (t_yaw_bLVLH_new(i-1)+k*360)-t_yaw_bLVLH_new(i)>170
            k = k-1;
        end
        t_yaw_bLVLH_new(i) = t_yaw_bLVLH_new(i) - k*360;        
        if t_roll_bLVLH_new(i)-(t_roll_bLVLH_new(i-1)+j*360)>170
            j = j+1;
        end
        if (t_roll_bLVLH_new(i-1)+j*360)-t_roll_bLVLH_new(i)>170
            j = j-1;
        end
        t_roll_bLVLH_new(i) = t_roll_bLVLH_new(i) -j*360;
    end
end

% Add noise due to sensor measurements
noise = 5/3600; % Pointing accuracy of the sensors [arcsecond]->[degree]
t_yaw_bLVLH_f = t_yaw_bLVLH_new(end)-rem(t_yaw_bLVLH_new(end),360) + noise;
t_pitch_bLVLH_f = t_pitch_bLVLH_new(end)-rem(t_pitch_bLVLH_new(end),360) + noise;
t_roll_bLVLH_f = t_roll_bLVLH_new(end)-rem(t_roll_bLVLH_new(end),360) + noise;
eul = eul2rotm([t_yaw_bLVLH_f,t_roll_bLVLH_f,t_pitch_bLVLH_f],'XYZ');
dt_q_bLVLH_i = dcm2quat(eul);

% Display
figure
plot3(t_reci_new(:,1),t_reci_new(:,2),t_reci_new(:,3))
title('Orbit')
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
grid on

figure
plot(t_t,t_w_new(:,1))
hold on
plot(t_t,t_w_new(:,2))
plot(t_t,t_w_new(:,3))
title('Angular Velocity of F_E_C_I to F_b')
xlabel('Time [s]')
ylabel('\omega [rad/s]')
legend('\omega_x','\omega_y','\omega_z')
grid on
hold off

figure
plot(t_t,t_eps_bLVLH_new(:,1))
hold on
plot(t_t,t_eps_bLVLH_new(:,2))
plot(t_t,t_eps_bLVLH_new(:,3))
plot(t_t,t_eta_bLVLH_new)
title('Quaternion of F_L_V_L_H to F_b')
xlabel('Time [s]')
ylabel('q')
grid on
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')
hold off
figure
plot(t_t,t_yaw_bLVLH_new)
hold on
plot(t_t,t_pitch_bLVLH_new)
plot(t_t,t_roll_bLVLH_new)
title('Euler Angles of F_L_V_L_H to F_b')
xlabel('Time [s]')
ylabel('Angles [^o]')
legend('\psi-yaw','\theta-pitch','\phi-roll')
grid on
hold off

%% De-tumble Phase
global Torques T_c

% LVLH to Body transformation
dt_C_bLVLH = quat2dcm(dt_q_bLVLH_i); % Quaternion to DCM

% ECI to LVLH transformation
dt_r = norm(t_reci); % Position in ECI
dt_z_hat = -t_reci/t_r; % Z LVLH
dt_y_hat = -cross(t_reci,t_veci)/norm(cross(t_reci,t_veci)); % Y LVLH
dt_x_hat = cross(dt_y_hat,dt_z_hat); % X LVLH
dt_C_LVLHG = [dt_x_hat,dt_y_hat,dt_z_hat]'; % ECI to LVLH DCM

% ECI to body transformation
dt_C_bG = dt_C_bLVLH*dt_C_LVLHG; % ECI to body DCM
dt_eta_bG = sqrt(trace(dt_C_bG)+1)/2; % Quaternion scalar

% Quaternion vectors
dt_eps_bG_1 = (dt_C_bG(2,3)-dt_C_bG(3,2))/(4*dt_eta_bG);
dt_eps_bG_2 = (dt_C_bG(3,1)-dt_C_bG(1,3))/(4*dt_eta_bG);
dt_eps_bG_3 = (dt_C_bG(1,2)-dt_C_bG(2,1))/(4*dt_eta_bG);
dt_eps_bG = [dt_eps_bG_1; dt_eps_bG_2; dt_eps_bG_3];

dt_zeta = 0.70; % Damping coefficient 70 [%]
dt_t_s = period;
dt_w_n = 4.4/(dt_t_s*dt_zeta); % Undamped natural frequency
%w_n = 0.5; % Natural frequency [rad/s]
dt_K_p = 2*(dt_w_n^2)*I_sc; % Proportional damping constant 
dt_K_d = 2*dt_zeta*dt_w_n*I_sc; % Derivative demping constant
q_c_bLVLH = [1;0;0;0]'; % Desired quaternion for LVLH to body
dt_eps_bLVLH_i = [dt_q_bLVLH_i(2);dt_q_bLVLH_i(3);dt_q_bLVLH_i(4)];
dt_eta_bLVLH_i = dt_q_bLVLH_i(1);
tinitial = 0; % Initial time [s]
tfinal = period*2; % Orbit time [s]
dt_t = [tinitial tfinal]; % Time span [s]
options = odeset('RelTol',1e-8,'AbsTol',1e-8,'OutputFcn',@outputfcndetumble); % Ode settings
dt_state = [t_reci;t_veci; dt_eps_bG; dt_eta_bG; dt_eps_bLVLH_i; dt_eta_bLVLH_i;t_w_bG_i]; % State vectors
[dt_t,dt_statenew]=ode45(@detumble,dt_t,dt_state,options,mu,I_sc,dt_K_p,dt_K_d,q_c_bLVLH);
% State vectors
dt_r_eci_new = [dt_statenew(:,1),dt_statenew(:,2),dt_statenew(:,3)];
dt_v_eci_new = [dt_statenew(:,4),dt_statenew(:,5),dt_statenew(:,6)];
% ECI to body
dt_eps_bG_new = [dt_statenew(:,7),dt_statenew(:,8),dt_statenew(:,9)];
dt_eta_bG_new = dt_statenew(:,10);
% LVLH to body
dt_eps_bLVLH_new = [dt_statenew(:,11),dt_statenew(:,12),dt_statenew(:,13)];
dt_eta_bLVLH_new = dt_statenew(:,14);
% Spacecraft angular velocity
dt_w_bG_sc_new = [dt_statenew(:,15),dt_statenew(:,16),dt_statenew(:,17)];
% Spacecraft angular momentum
for i = 1:length(dt_t)
    dt_h_sc_new_1(i,:) = I_sc*dt_w_bG_sc_new(i,:)';
    dt_h_sc_new_mag_1(i) = norm(dt_h_sc_new_1(i,:));
end
% ECI to body Euler angles
k = 0;
j = 0;
for i = 1:length(dt_t)
    dt_C_bG_new = quat2dcm([dt_eta_bG_new(i),dt_eps_bG_new(i,:)]);
    [dt_yaw_bG_new(i),dt_pitch_bG_new(i),dt_roll_bG_new(i)] = dcm_to_ypr(dt_C_bG_new);

    %Angle fixing beyond 180 due matlab atan2d restrictions
    if i>1
        if dt_yaw_bG_new(i)-(dt_yaw_bG_new(i-1)+k*360)>170
            k = k+1;
        end
        if (dt_yaw_bG_new(i-1)+k*360)-dt_yaw_bG_new(i)>170
            k = k-1;
        end
        dt_yaw_bG_new(i) = dt_yaw_bG_new(i) - k*360;        
        if dt_roll_bG_new(i)-(dt_roll_bG_new(i-1)+j*360)>170
            j = j+1;
        end
        if (dt_roll_bG_new(i-1)+j*360)-dt_roll_bG_new(i)>170
            j = j-1;
        end
        dt_roll_bG_new(i) = dt_roll_bG_new(i) -j*360;
    end
end
k = 0;
j = 0;
for i = 1:length(dt_t)
    dt_C_bLVLH_new = quat2dcm([dt_eta_bLVLH_new(i),dt_eps_bLVLH_new(i,:)]);
    [dt_yaw_bLVLH_new(i),dt_pitch_bLVLH_new(i),dt_roll_bLVLH_new(i)] = dcm_to_ypr(dt_C_bLVLH_new);

    %Angle fixing beyond 180 due matlab atan2d restrictions
    if i>1
        if dt_yaw_bLVLH_new(i)-(dt_yaw_bLVLH_new(i-1)+k*360)>170
            k = k+1;
        end
        if (dt_yaw_bLVLH_new(i-1)+k*360)-dt_yaw_bLVLH_new(i)>170
            k = k-1;
        end
        dt_yaw_bLVLH_new(i) = dt_yaw_bLVLH_new(i) - k*360;        
        if dt_roll_bLVLH_new(i)-(dt_roll_bLVLH_new(i-1)+j*360)>170
            j = j+1;
        end
        if (dt_roll_bLVLH_new(i-1)+j*360)-dt_roll_bLVLH_new(i)>170
            j = j-1;
        end
        dt_roll_bLVLH_new(i) = dt_roll_bLVLH_new(i) -j*360;
    end
end

% Accumulated Torque
for i = 1:length(Torques(:,1))
    if i == 1
    accu_T_x(i) = abs(Torques(i,1));
    accu_T_y(i) = abs(Torques(i,2));
    accu_T_z(i) = abs(Torques(i,3));
    else
    accu_T_x(i) = abs(Torques(i,1)) + accu_T_x(i-1);
    accu_T_y(i) = abs(Torques(i,2)) + accu_T_x(i-1);
    accu_T_z(i) = abs(Torques(i,3)) + accu_T_x(i-1);
    end
end
g = 9.81;
Isp = 3200;
d = 0.75*2;
mprop_x = accu_T_x(end)*dt_t_s/(g*Isp*d);
mprop_y = accu_T_y(end)*dt_t_s/(g*Isp*d);
mprop_z = accu_T_z(end)*dt_t_s/(g*Isp*d);
mprop_tot = (mprop_x + mprop_y + mprop_z)*10^3; % micrograms

fprintf('Total fuel mass required: %.3f [gram]',mprop_tot)

figure
plot(dt_t,dt_w_bG_sc_new(:,1))
hold on
plot(dt_t,dt_w_bG_sc_new(:,2))
plot(dt_t,dt_w_bG_sc_new(:,3))
title('Angular Velocity of Spacecraft')
xlabel('Time [s]')
ylabel('\omega_S_/_C [rad/s]')
legend('\omega_x','\omega_y','\omega_z')
grid on
hold off

figure
plot(dt_t,dt_h_sc_new_mag_1)
title('Total Angular Momentum Accumulated of Spacecraft')
xlabel('Time [s]')
ylabel('Angular Momentum [N*m/s]')
grid on


figure
plot(dt_t,dt_eps_bLVLH_new(:,1))
hold on
plot(dt_t,dt_eps_bLVLH_new(:,2))
plot(dt_t,dt_eps_bLVLH_new(:,3))
plot(dt_t,dt_eta_bLVLH_new(:,1))
title('Quaternion of F_L_V_L_H to F_b')
xlabel('Time [s]')
ylabel('q')
grid on
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')
hold off

figure
plot(dt_t,dt_yaw_bLVLH_new)
hold on
plot(dt_t,dt_pitch_bLVLH_new)
plot(dt_t,dt_roll_bLVLH_new)
title('Euler Angles of F_L_V_L_H to F_b')
xlabel('Time [s]')
ylabel('Angles [^o]')
legend('\psi-yaw','\theta-pitch','\phi-roll')
grid on
hold off

figure
plot(dt_t,dt_eps_bG_new(:,1))
hold on
plot(dt_t,dt_eps_bG_new(:,2))
plot(dt_t,dt_eps_bG_new(:,3))
plot(dt_t,dt_eta_bG_new(:,1))
title('Quaternion of F_E_C_I to F_b')
xlabel('Time [s]')
ylabel('q')
grid on
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')
hold off

figure
plot(dt_t,dt_yaw_bG_new)
hold on
plot(dt_t,dt_pitch_bG_new)
plot(dt_t,dt_roll_bG_new)
title('Euler Angles of F_E_C_I to F_b')
xlabel('Time [s]')
ylabel('Angles [^o]')
legend('\psi-yaw','\theta-pitch','\phi-roll')
grid on
hold off

figure
plot(Torques(:,4),Torques(:,1))
hold on
plot(Torques(:,4),Torques(:,2))
plot(Torques(:,4),Torques(:,3))
title('Control Torque')
xlabel('Time [s]')
ylabel('Torque [N*m]')
legend('T_x','T_y','T_z')
grid on
hold off

figure
plot(Torques(:,4),accu_T_x)
hold on
plot(Torques(:,4),accu_T_y)
plot(Torques(:,4),accu_T_z)
title('Accumulated Control Torque')
xlabel('Time [s]')
ylabel('Torque [N*m]')
legend('T_x','T_y','T_z')
grid on
hold off

% 
% % % Animation for fun
% figure
% title('Spacecraft Motion Simulation')
% hold on
% [M] = RBMotionSim(dt_yaw_bLVLH_new, dt_pitch_bLVLH_new, dt_roll_bLVLH_new);
% hold off

%% Operational Phase
global Torques T T_ad T_sp T_gg T_m T_c

%% Moment of inertia calculation post deployment
% Center of mass of components
cm_solar_1 = [0;-2.5;0]; % Solar panel 1 [m]
cm_solar_2 = [0;2.5;0]; % Solar panel 2 [m]
cm_sensor = [0;0;1.5]; % Sensor [m]
cm_cube = [0;0;0]; % Bus [m]
% Center of mass of spacecraft [m]
cm_sc = (m_solar*cm_solar_1+m_solar*cm_solar_2+m_sensor*cm_sensor+m_cube*cm_cube)/m_sc_tot;

% Distance to center of mass of spacecraft from each component
dist_solar_1 = cm_solar_1-cm_sc; % Solar panel 1 [m]
dist_solar_2 = cm_solar_2-cm_sc; % Solar panel 2 [m]
dist_sensor = cm_sensor-cm_sc; % Sensor [m]
dist_cube = cm_cube-cm_sc; % Bus [m]

% Crossed distances to use for parallel axis theorem
dist_solar_1_cross = cross_x(dist_solar_1(1),dist_solar_1(2),dist_solar_1(3)); % Solar panel 1 [m]
dist_solar_2_cross = cross_x(dist_solar_2(1),dist_solar_2(2),dist_solar_2(3)); % Solar panel 2 [m]
dist_sensor_cross = cross_x(dist_sensor(1),dist_sensor(2),dist_sensor(3)); % Sensor [m]
dist_cube_cross = cross_x(dist_cube(1),dist_cube(2),dist_cube(3)); % Bus [m]

% Moment of inertia function for rectangular shapes
I_rect = @(m,y,x,z) m/12*[x^2+z^2 0 0 ; 0 y^2+z^2 0; 0 0 y^2+x^2];

% Moment of inertia for each component
I_solar_1 = I_rect(m_solar,h_solar,w_solar,d_solar) - dist_solar_1_cross*dist_solar_1_cross*m_solar; % Solar panel 1 [kg*m^2]
I_solar_2 = I_rect(m_solar,h_solar,w_solar,d_solar) - dist_solar_2_cross*dist_solar_2_cross*m_solar; % Solar panel 2 [kg*m^2]
I_sensor = I_rect(m_sensor,h_sensor,w_sensor,d_sensor) - dist_sensor_cross*dist_sensor_cross*m_sensor; % Sensor [kg*m^2]
I_cube = I_sc - dist_cube_cross*dist_cube_cross*m_cube; % Bus [kg*m^2]

% Moment of inertia for spacecraft [kg*m^2]
I_sc_tot = I_solar_1 + I_solar_2 + I_sensor +I_cube;

% 3-axis stabilize Control law
no_zeta = 0.70; % Damping coefficient 70 [%]
no_t_s = 100; % Settling time [s]
no_w_n = 4.4/(no_t_s*no_zeta); % Undamped natural frequency [rad/s]
%w_n = 0.5; % Natural frequency 
no_K_p = 2*(no_w_n^2)*I_sc_tot; % Proportional damping constant 
no_K_d = 2*no_zeta*no_w_n*I_sc_tot; % Derivative demping constant
% Wheel moment of inertia
I_ws = 0.037; % Spin axis [kg*m^2]
I_wt = 0.0769; % Transverse axis [kg*m^2]
m_w = 4.85; % Mass of wheel [kg]

% Spacecraft moment of inertia with wheels [kg*m^2]
I_sc_w = (I_sc_tot + (2*I_wt+I_ws+2*m_w)*eye(3));
% Wheel moment of inetia [kg*m^2]
I_w = [I_ws 0 0;
    0 I_ws 0;
    0 0 I_ws];
q_c_bLVLH = [1;0;0;0]'; % Desired quaternion for LVLH to body
no_w_bG_w_i = [0;0;0]; % Initial wheel angular velocity [rad/s]
no_tfinal = period*5; % Orbit time [s]
no_t = [tinitial no_tfinal]; % Time span [s]
options = odeset('RelTol',1e-8,'AbsTol',1e-8,'OutputFcn',@outputfcnFinal); % Ode settings
no_state = [t_reci;t_veci; t_eps_bG; t_eta_bG; t_eps_bLVLH_i; t_eta_bLVLH_i;no_w_bG_sc_i;no_w_bG_w_i]; % State vectors
[no_t,no_statenew]=ode45(@normal_operation,no_t,no_state,options,mu,a_x,a_y,a_z,a_solar,a_cube,a_sensor_side,a_sensor_bottom,sp,sun_ECI,C_D, mdp,I_sc_w,I_w,no_K_p,no_K_d,q_c_bLVLH);
% State vectors
no_r_eci_new = [no_statenew(:,1),no_statenew(:,2),no_statenew(:,3)];
no_v_eci_new = [no_statenew(:,4),no_statenew(:,5),no_statenew(:,6)];
% ECI to body
no_eps_bG_new = [no_statenew(:,7),no_statenew(:,8),no_statenew(:,9)];
no_eta_bG_new = no_statenew(:,10);
% LVLH to body
no_eps_bLVLH_new = [no_statenew(:,11),no_statenew(:,12),no_statenew(:,13)];
no_eta_bLVLH_new = no_statenew(:,14);
% Spacecraft angular velocity
no_w_bG_sc_new = [no_statenew(:,15),no_statenew(:,16),no_statenew(:,17)];
% Wheel angulr velocity
no_w_bG_w_new = [no_statenew(:,18),no_statenew(:,19),no_statenew(:,20)];
% Spacecraft angular momentum
for i = 1:length(no_t)
    if i == 1
    no_h_sc_new_1(i,:) = I_sc_tot*no_w_bG_sc_new(i,:)';
    no_h_sc_new_x(i) = no_h_sc_new_1(i,1);
    no_h_sc_new_y(i) = no_h_sc_new_1(i,2);
    no_h_sc_new_z(i) = no_h_sc_new_1(i,3);
    no_h_sc_new_mag_1(i) = norm(no_h_sc_new_1(i,:));
    else
    no_h_sc_new_1(i,:) = I_sc_tot*no_w_bG_sc_new(i,:)';
    no_h_sc_new_x(i) = no_h_sc_new_1(i,1)+no_h_sc_new_x(i-1);
    no_h_sc_new_y(i) = no_h_sc_new_1(i,2)+no_h_sc_new_y(i-1);
    no_h_sc_new_z(i) = no_h_sc_new_1(i,3)+no_h_sc_new_z(i-1);
    no_h_sc_new_mag_1(i) = norm(no_h_sc_new_1(i,:));
    end
end
% Wheel angular momentum
for i = 1:length(no_t)
    if i == 1
    no_h_w_new_1(:,i) = I_w*no_w_bG_w_new(i,:)';
    no_h_w_new_x(i) = no_h_w_new_1(1,i);
    no_h_w_new_y(i) = no_h_w_new_1(2,i);
    no_h_w_new_z(i) = no_h_w_new_1(3,i);
    no_h_w_new_mag_1(i) = norm(no_h_w_new_1(:,i));
    else
    no_h_w_new_1(:,i) = I_w*no_w_bG_w_new(i,:)';
    no_h_w_new_x(i) = no_h_w_new_1(1,i)+no_h_w_new_x(i-1);
    no_h_w_new_y(i) = no_h_w_new_1(2,i)+no_h_w_new_y(i-1);
    no_h_w_new_z(i) = no_h_w_new_1(3,i)+no_h_w_new_z(i-1);
    no_h_w_new_mag_1(i) = norm(no_h_w_new_1(:,i));
    end
end


% ECI to body Euler angles
k = 0;
j = 0;
for i = 1:length(no_t)
    no_C_bG_new = quat2dcm([no_eta_bG_new(i),no_eps_bG_new(i,:)]);
    [no_yaw_bG_new(i),no_pitch_bG_new(i),no_roll_bG_new(i)] = dcm_to_ypr(no_C_bG_new);

    %Angle fixing beyond 180 due matlab atan2d restrictions
    if i>1
        if no_yaw_bG_new(i)-(no_yaw_bG_new(i-1)+k*360)>170
            k = k+1;
        end
        if (no_yaw_bG_new(i-1)+k*360)-no_yaw_bG_new(i)>170
            k = k-1;
        end
        no_yaw_bG_new(i) = no_yaw_bG_new(i) - k*360;        
        if no_roll_bG_new(i)-(no_roll_bG_new(i-1)+j*360)>170
            j = j+1;
        end
        if (no_roll_bG_new(i-1)+j*360)-no_roll_bG_new(i)>170
            j = j-1;
        end
        no_roll_bG_new(i) = no_roll_bG_new(i) -j*360;
    end
end
k = 0;
j = 0;
for i = 1:length(no_t)
    no_C_bLVLH_new = quat2dcm([no_eta_bLVLH_new(i),no_eps_bLVLH_new(i,:)]);
    [no_yaw_bLVLH_new(i),no_pitch_bLVLH_new(i),no_roll_bLVLH_new(i)] = dcm_to_ypr(no_C_bLVLH_new);

    %Angle fixing beyond 180 due matlab atan2d restrictions
    if i>1
        if no_yaw_bLVLH_new(i)-(no_yaw_bLVLH_new(i-1)+k*360)>170
            k = k+1;
        end
        if (no_yaw_bLVLH_new(i-1)+k*360)-no_yaw_bLVLH_new(i)>170
            k = k-1;
        end
        no_yaw_bLVLH_new(i) = no_yaw_bLVLH_new(i) - k*360;        
        if no_roll_bLVLH_new(i)-(no_roll_bLVLH_new(i-1)+j*360)>170
            j = j+1;
        end
        if (no_roll_bLVLH_new(i-1)+j*360)-no_roll_bLVLH_new(i)>170
            j = j-1;
        end
        no_roll_bLVLH_new(i) = no_roll_bLVLH_new(i) -j*360;
    end
end

figure
plot(no_t,no_w_bG_sc_new(:,1))
hold on
plot(no_t,no_w_bG_sc_new(:,2))
plot(no_t,no_w_bG_sc_new(:,3))
title('Angular Velocity of Spacecraft')
xlabel('Time [s]')
ylabel('\omega_S_/_C [rad/s]')
legend('\omega_x','\omega_y','\omega_z')
grid on
hold off

figure
plot(no_t,no_h_sc_new_1(:,1))
hold on
plot(no_t,no_h_sc_new_1(:,2))
plot(no_t,no_h_sc_new_1(:,3))
plot(no_t,no_h_sc_new_mag_1)
title('Angular Momentum of Spacecraft')
xlabel('Time [s]')
ylabel('Angular Momentum [N*m/s]')
legend('x','y','z','norm')
grid on
hold off

figure
plot(no_t,no_h_sc_new_x)
hold on
plot(no_t,no_h_sc_new_y)
plot(no_t,no_h_sc_new_z)
title('Total Angular Momentum Accumulated of Spacecraft')
xlabel('Time [s]')
ylabel('Angular Momentum [N*m/s]')
legend('x','y','z')
grid on
hold off

figure
plot(no_t,no_w_bG_w_new(:,1))
hold on
plot(no_t,no_w_bG_w_new(:,2))
plot(no_t,no_w_bG_w_new(:,3))
title('Angular Velocity of Wheel')
xlabel('Time [s]')
ylabel('\omega_W [rad/s]')
legend('\omega_x','\omega_y','\omega_z')
grid on
hold off

figure
plot(no_t,no_h_w_new_1(1,:))
hold on
plot(no_t,no_h_w_new_1(2,:))
plot(no_t,no_h_w_new_1(3,:))
plot(no_t,no_h_w_new_mag_1)
title('Angular Momentum of Wheel')
xlabel('Time [s]')
ylabel('Angular Momentum [N*m/s]')
legend('x','y','z','norm')
grid on
hold off

figure
plot(no_t,no_h_w_new_x)
hold on
plot(no_t,no_h_w_new_y)
plot(no_t,no_h_w_new_z)
title('Total Angular Momentum Accumulated of Wheel')
xlabel('Time [s]')
ylabel('Angular Momentum [N*m/s]')
legend('x','y','z')
grid on
hold off

figure
plot(no_t,no_eps_bLVLH_new(:,1))
hold on
plot(no_t,no_eps_bLVLH_new(:,2))
plot(no_t,no_eps_bLVLH_new(:,3))
plot(no_t,no_eta_bLVLH_new(:,1))
title('Quaternion of F_L_V_L_H to F_b')
xlabel('Time [s]')
ylabel('q')
grid on
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')
hold off

figure
plot(no_t,no_yaw_bLVLH_new)
hold on
plot(no_t,no_pitch_bLVLH_new)
plot(no_t,no_roll_bLVLH_new)
title('Euler Angles of F_L_V_L_H to F_b')
xlabel('Time [s]')
ylabel('Angles [^o]')
legend('\psi-yaw','\theta-pitch','\phi-roll')
grid on
hold off

figure
plot(no_t,no_eps_bG_new(:,1))
hold on
plot(no_t,no_eps_bG_new(:,2))
plot(no_t,no_eps_bG_new(:,3))
plot(no_t,no_eta_bG_new(:,1))
title('Quaternion of F_E_C_I to F_b')
xlabel('Time [s]')
ylabel('q')
grid on
legend('\epsilon_1','\epsilon_2','\epsilon_3','\eta')
hold off

figure
plot(no_t,no_yaw_bG_new)
hold on
plot(no_t,no_pitch_bG_new)
plot(no_t,no_roll_bG_new)
title('Euler Angles of F_E_C_I to F_b')
xlabel('Time [s]')
ylabel('Angles [^o]')
legend('\psi-yaw','\theta-pitch','\phi-roll')
grid on
hold off

figure
subplot(2,3,1)
plot(Torques(:,19),Torques(:,1))
hold on
plot(Torques(:,19),Torques(:,2))
plot(Torques(:,19),Torques(:,3))
title('Total Torque')
xlabel('Time [s]')
ylabel('Torque [N*m]')
legend('T_x','T_y','T_z')
grid on
hold off
subplot(2,3,2)
plot(Torques(:,19),Torques(:,4))
hold on
plot(Torques(:,19),Torques(:,5))
plot(Torques(:,19),Torques(:,6))
title('Aerodynamic Torque NRLMSISE00')
xlabel('Time [s]')
ylabel('Torque [N*m]')
legend('T_x','T_y','T_z')
grid on
hold off
subplot(2,3,3)
plot(Torques(:,19),Torques(:,7))
hold on
plot(Torques(:,19),Torques(:,8))
plot(Torques(:,19),Torques(:,9))
title('Solar Pressure Torque')
xlabel('Time [s]')
ylabel('Torque [N*m]')
legend('T_x','T_y','T_z')
grid on
hold off
subplot(2,3,4)
plot(Torques(:,19),Torques(:,10))
hold on
plot(Torques(:,19),Torques(:,11))
plot(Torques(:,19),Torques(:,12))
title('Gravity Gradient Torque')
xlabel('Time [s]')
ylabel('Torque [N*m]')
legend('T_x','T_y','T_z')
grid on
hold off
subplot(2,3,5)
plot(Torques(:,19),Torques(:,13))
hold on
plot(Torques(:,19),Torques(:,14))
plot(Torques(:,19),Torques(:,15))
title('Magnetic Torque')
xlabel('Time [s]')
ylabel('Torque [N*m]')
legend('T_x','T_y','T_z')
grid on
hold off
subplot(2,3,6)
plot(Torques(:,19),Torques(:,16))
hold on
plot(Torques(:,19),Torques(:,17))
plot(Torques(:,19),Torques(:,18))
title('Control Torque')
xlabel('Time [s]')
ylabel('Torque [N*m]')
legend('T_x','T_y','T_z')
grid on
hold off
sgtitle('All Torques')

%% Functions
%% COE to RV
% Ryo Takatori

function [ rperi, vperi, reci, veci ] = coe2rv( mu, e, h, inc, OMEGA, omega, theta)
% State vectors in perifocal
rperi = ((h^2)/mu)*(1/(1+e*cosd(theta)))*[cosd(theta) sind(theta) 0];
vperi = (mu/h)*[-sind(theta) e+cosd(theta) 0];

% Perifocal to ECI
[ reci, veci ] = peri2eci( rperi, vperi, omega, OMEGA, inc );

end

%% Perifocal to ECI
function [ reci, veci ] = peri2eci( rperi, vperi, omega, OMEGA, inc )
% Rotate to ECI
R3omega = [cosd(omega) sind(omega) 0;
    -sind(omega) cosd(omega) 0;
    0 0 1];

R1inc = [1 0 0;
    0 cosd(inc) sind(inc);
    0 -sind(inc) cosd(inc)];

R3OMEGA = [cosd(OMEGA) sind(OMEGA) 0;
    -sind(OMEGA) cosd(OMEGA) 0;
    0 0 1];

Q = (R3omega*R1inc*R3OMEGA)';

reci = Q*rperi';
veci = Q*vperi';
end

%% Tumble
function [dstatedt] = tumble(t,state,mu,I)
%% Two body
dx = state(4); % dx [km/s]
dy = state(5); % dy [km/s]
dz = state(6); % dz [km/s]

r = norm([state(1) state(2) state(3)]); % Magnitude of position [km]
reci = [state(1);state(2);state(3)]; % Positon vector in ECI [km]
veci = [state(4); state(5); state(6)]; % Velocity vector in ECI [km/s]
% Equation of motion
ddx = -mu*state(1)/r^3; % ddx [km/s^2]
ddy = -mu*state(2)/r^3; % ddy [km/s^2]
ddz = -mu*state(3)/r^3; % ddz [km/s^2]
% ECI to body
eps_bG = [state(7); state(8); state(9)]; % Quaternion vector
eta_bG = state(10); % Quaternion scalar

% LVLH to body
eps_bLVLH = [state(11); state(12); state(13)];
eta_bLVLH = state(14);
C_bLVLH = quat2dcm([eta_bLVLH;eps_bLVLH]');

% ECI to LVLH
z_hat = -reci/r; % Z LVLH
y_hat = -cross(reci,veci)/norm(cross(reci,veci)); % Y LVLH
x_hat = cross(y_hat,z_hat); % X LVLH
C_LVLHG = [x_hat,y_hat,z_hat]'; % ECI to LVLH DCM 

% ECI to body matrix
C_bG = C_bLVLH*C_LVLHG; % ECI to body DCM

% Angular velocities in various frames
omega_bG = [state(15);state(16);state(17)]; % Angular velocity [rad/s]
omega_LVLHG = C_bG*(cross(reci,veci)/r^2);
omega_bLVLH = omega_bG-omega_LVLHG;

q_bG = [eta_bG;eps_bG]; % Quaternion in MATLAB form
q_bLVLH = [eta_bLVLH;eps_bLVLH]; % Quaternion in MATLAB form
    
deps_bG = 0.5*(q_bG(1)*eye(3)+cross_x(q_bG(2),q_bG(3),q_bG(4)))*omega_bG; % Derivative of quaternion vector
deta_bG = -0.5*[q_bG(2);q_bG(3);q_bG(4)]'*omega_bG; % Derivative of quaternion scalar
deps_bLVLH = 0.5*(q_bLVLH(1)*eye(3)+cross_x(q_bLVLH(2),q_bLVLH(3),q_bLVLH(4)))*omega_bLVLH; % Quaternion vector derivative
deta_bLVLH = -0.5*[q_bLVLH(2);q_bLVLH(3);q_bLVLH(4)]'*omega_bLVLH; % Quaternion scalar derivative
domega_bG = (I^-1)*(-cross_x(omega_bG(1),omega_bG(2),omega_bG(3))*I*omega_bG); % Derivative of angular velocity [rad/s^2]
dstatedt = [dx;dy;dz;ddx;ddy;ddz;deps_bG;deta_bG;deps_bLVLH;deta_bLVLH;domega_bG];

end

%% DCM to YPR
function [yaw pitch roll] = dcm_to_ypr(Q)
%{
This function finds the angles of the yaw-pitch-roll sequence
R1(gamma)*R2(beta)*R3(alpha) from the direction cosine matrix.
Q - direction cosine matrix
yaw - yaw angle (deg)
pitch - pitch angle (deg)
roll - roll angle (deg)
User M-function required: atan2d_0_360
%}
yaw = atan2d(Q(1,2), Q(1,1));
pitch = -asind(Q(1,3));
roll = atan2d(Q(2,3), Q(3,3));
end

%% De-Tumble
function [dstatedt] = detumble(t, state,mu,I_sc,K_p,K_d,q_c_bLVLH)

%% Two body
dx = state(4); % dx [km/s]
dy = state(5); % dy [km/s]
dz = state(6); % dz [km/s]

r = norm([state(1) state(2) state(3)]); % Magnitude of position [km]
reci = [state(1);state(2);state(3)]; % Positon vector in ECI [km]
veci = [state(4); state(5); state(6)]; % Velocity vector in ECI [km/s]
% Equation of motion
ddx = -mu*state(1)/r^3; % ddx [km/s^2]
ddy = -mu*state(2)/r^3; % ddy [km/s^2]
ddz = -mu*state(3)/r^3; % ddz [km/s^2]

%% Quaternion and Rotation Matrices
% ECI to body
eps_bG = [state(7); state(8); state(9)]; % Quaternion vector
eta_bG = state(10); % Quaternion scalar

% LVLH to body
eps_bLVLH = [state(11); state(12); state(13)];
eta_bLVLH = state(14);
C_bLVLH = quat2dcm([eta_bLVLH;eps_bLVLH]');

% ECI to LVLH
z_hat = -reci/r; % Z LVLH
y_hat = -cross(reci,veci)/norm(cross(reci,veci)); % Y LVLH
x_hat = cross(y_hat,z_hat); % X LVLH
C_LVLHG = [x_hat,y_hat,z_hat]'; % ECI to LVLH DCM 

% ECI to body matrix
C_bG = C_bLVLH*C_LVLHG; % ECI to body DCM

% Angular velocities in various frames
omega_bG = [state(15);state(16);state(17)]; % Angular velocity [rad/s]
omega_LVLHG = C_bG*(cross(reci,veci)/r^2);
omega_bLVLH = omega_bG-omega_LVLHG;


 % Quaternion errors
    q_bG = [eta_bG;eps_bG]; % Quaternion in MATLAB form
    q_bLVLH = [eta_bLVLH;eps_bLVLH]; % Quaternion in MATLAB form
    q_e_bLVLH = quatmultiply(quatconj(q_c_bLVLH),q_bLVLH'); % Quaternion error
    T_c = -K_p*[q_e_bLVLH(2);q_e_bLVLH(3);q_e_bLVLH(4)]-K_d*omega_bLVLH;
    deps_bG = 0.5*(q_bG(1)*eye(3)+cross_x(q_bG(2),q_bG(3),q_bG(4)))*omega_bG; % Derivative of quaternion vector
    deta_bG = -0.5*[q_bG(2);q_bG(3);q_bG(4)]'*omega_bG; % Derivative of quaternion scalar
    deps_bLVLH = 0.5*(q_e_bLVLH(1)*eye(3)+cross_x(q_e_bLVLH(2),q_e_bLVLH(3),q_e_bLVLH(4)))*omega_bLVLH; % Quaternion vector derivative
    deta_bLVLH = -0.5*[q_e_bLVLH(2);q_e_bLVLH(3);q_e_bLVLH(4)]'*omega_bLVLH; % Quaternion scalar derivative
    domega_bG = (I_sc^-1)*(-cross_x(omega_bG(1),omega_bG(2),omega_bG(3))*I_sc*omega_bG+T_c); % Derivative of angular velocity [rad/s^2]
    dstatedt = [dx;dy;dz;ddx;ddy;ddz;deps_bG;deta_bG;deps_bLVLH;deta_bLVLH;domega_bG]; % Return state vectors
end

function [dstatedt] = normal_operation(t, state, mu,a_x,a_y,a_z,a_solar,a_cube,a_sensor_side,a_sensor_bottom,sp,sun,C_D, mdp,I_sc,I_w,K_p,K_d,q_c_bLVLH)

%% Two body
dx = state(4); % dx [km/s]
dy = state(5); % dy [km/s]
dz = state(6); % dz [km/s]

r = norm([state(1) state(2) state(3)]); % Magnitude of position [km]
reci = [state(1);state(2);state(3)]; % Positon vector in ECI [km]
veci = [state(4); state(5); state(6)]; % Velocity vector in ECI [km/s]
% Equation of motion
ddx = -mu*state(1)/r^3; % ddx [km/s^2]
ddy = -mu*state(2)/r^3; % ddy km/s^2]
ddz = -mu*state(3)/r^3; % ddz [km/s^2]

%% Quaternion and Rotation Matrices
% ECI to body
eps_bG = [state(7); state(8); state(9)]; % Quaternion vector
eta_bG = state(10); % Quaternion scalar

% LVLH to body
eps_bLVLH = [state(11); state(12); state(13)];
eta_bLVLH = state(14);
C_bLVLH = quat2dcm([eta_bLVLH;eps_bLVLH]');

% ECI to LVLH
z_hat = -reci/r; % Z LVLH
y_hat = -cross(reci,veci)/norm(cross(reci,veci)); % Y LVLH
x_hat = cross(y_hat,z_hat); % X LVLH
C_LVLHG = [x_hat,y_hat,z_hat]'; % ECI to LVLH DCM 

% ECI to body matrix
C_bG = C_bLVLH*C_LVLHG; % ECI to body DCM

% Angular velocities in various frames
omega_bG_sc = [state(15);state(16);state(17)]; % Angular velocity [rad/s]
omega_LVLHG = C_bG*(cross(reci,veci)/r^2);
omega_bLVLH = omega_bG_sc-omega_LVLHG;
omega_bG_w = [state(18);state(19);state(20)];

% Position vectors
r_b = C_bG*reci; % Position in body [km]
v_b = C_bG*veci; % Velocity in body [km/s]

%% Spacecraft Properties
% Total Spacecraft
n_x_b_unit = [1;0;0]; % n vector for x face
n_y_b_unit = [0;1;0]; % n vector for y face
n_z_b_unit = [0;0;1]; % n vector for z face
n_x_b_unit_neg = [-1;0;0]; % n vector for x face
n_y_b_unit_neg = [0;-1;0]; % n vector for y face
n_z_b_unit_neg = [0;0;-1]; % n vector for z face

%% Solar pressure
sun_b = (C_bG*sun)/(norm(C_bG*sun)); % Unit sun vector in body
n_x_s = dot(n_x_b_unit,sun_b); % amount of x area
n_y_s = dot(n_y_b_unit,sun_b); % amount of y area
n_z_s = dot(n_z_b_unit,sun_b); % amount of z
n_x_s_neg = dot(n_x_b_unit_neg,sun_b); % amount of x area
n_y_s_neg = dot(n_y_b_unit_neg,sun_b); % amount of y area
n_z_s_neg = dot(n_z_b_unit_neg,sun_b); % amount of z area
% Checking for direction and applicable area
% Rho vectos from the center of mass to center of surface
% Cube
rho_cube_x_s = [1;0;0]*a_cube;
rho_cube_y_s = [0;1;0]*a_cube;
rho_cube_z_s = [0;0;1]*a_cube;
rho_cube_x_s_neg = [-1;0;0]*a_cube;
rho_cube_y_s_neg = [0;-1;0]*a_cube;
rho_cube_z_s_neg = [0;0;-1]*a_cube;

% Solar panel 1
rho_solar_z_s_1 = [0;2.5;0]*a_solar;
rho_solar_z_s_neg_1 = [0;2.5;0]*a_solar;

% Solar panel 2
rho_solar_z_s_2 = [0;-2.5;0]*a_solar;
rho_solar_z_s_neg_2 = [0;-2.5;0]*a_solar;

% Sensor
rho_sensor_x_s = [0.125;0;2]*a_sensor_side;
rho_sensor_y_s = [0;0.125;2]*a_sensor_side;
rho_sensor_z_s = [0;0;3]*a_sensor_bottom;
rho_sensor_x_s_neg = [-0.125;0;2]*a_sensor_side;
rho_sensor_y_s_neg = [0;-0.125;2]*a_sensor_side;

% Checking for direction and applicable area
if n_x_s<0
    n_x_s = 0;
    rho_cube_x_s = [0;0;0];
    rho_sensor_x_s = [0;0;0];
end
if n_y_s<0
    n_y_s = 0;
    rho_cube_y_s = [0;0;0];
    rho_sensor_y_s = [0;0;0];
end
if n_z_s<0
    n_z_s = 0;
    rho_cube_z_s = [0;0;0];
    rho_sensor_z_s = [0;0;0];
    rho_solar_z_s_1 = [0;0;0];
    rho_solar_z_s_2 = [0;0;0];

end
if n_x_s_neg<0
    n_x_s_neg = 0;
    rho_cube_x_s_neg = [0;0;0];
    rho_sensor_x_s_neg = [0;0;0];
end
if n_y_s_neg<0
    n_y_s_neg = 0;
    rho_cube_y_s_neg = [0;0;0];
    rho_sensor_y_s_neg = [0;0;0];
end
if n_z_s_neg<0
    n_z_s_neg = 0;
    rho_cube_z_s_neg = [0;0;0];
    rho_solar_z_s_neg_1 = [0;0;0];
    rho_solar_z_s_neg_2 = [0;0;0];
end
rho_tot_s = (rho_cube_x_s*n_x_s + ...
    rho_cube_y_s*n_y_s + ...
    rho_cube_z_s*n_z_s + ...
    rho_cube_x_s_neg*n_x_s_neg + ...
    rho_cube_y_s_neg*n_y_s_neg + ...
    rho_cube_z_s_neg*n_z_s_neg + ...
    rho_solar_z_s_1*n_z_s + ...
    rho_solar_z_s_neg_1*n_z_s_neg + ...
    rho_solar_z_s_2*n_z_s + ...
    rho_solar_z_s_neg_2*n_z_s_neg + ...
    rho_sensor_x_s*n_x_s + ...
    rho_sensor_y_s*n_y_s + ...
    rho_sensor_z_s*n_z_s + ...
    rho_sensor_x_s_neg*n_x_s_neg + ...
    rho_sensor_y_s_neg*n_y_s_neg);
a_wet = a_x*(n_x_s+n_x_s_neg)+a_y*(n_y_s+n_y_s_neg)+a_z*(n_z_s+n_z_s_neg); % Wetted area [m^2]
c_sp = (rho_tot_s/a_wet); % Center of solar pressure
F_s = -(sp*a_wet*sun_b); % Force due to sun [kN]
T_sp = cross(c_sp,F_s); % Solar pressure torque [kN*km]

%% Gravity gradient
T_gg = ((3*mu)/norm(r_b)^5)*cross_x(r_b(1),r_b(2),r_b(3))*I_sc*r_b; % Gravity gradient torque [N*km]

%% Magnetic torque
% Time adjustments
ht = t/3600;
hrs = floor(ht);
mt = (ht-floor(ht))*60;
mins = floor(mt);
sec = (mt-floor(mt))*60;

decy = decyear(2019,5,31,12+hrs,mins,sec); % decimal year
m_b = mdp; % Magnetic dipole [A*km^2]
rlla = eci2lla(reci'*1000,[2019,5,31,12+hrs,mins,sec]); % Position in LLA [km]
[b_ned, h, dec, dip, f] = wrldmagm(rlla(3), rlla(1), rlla(2), decy, '2015'); % Magnetic field vector in NED [nT]
[b_ecef_x,b_ecef_y,b_ecef_z] = ned2ecef(b_ned(1),b_ned(2),b_ned(3),rlla(1),rlla(2),rlla(3),referenceSphere); % Magnetic field vector in ECEF [nT]
b_eci = dcmeci2ecef('IAU-2000/2006',[2019,5,31,12+hrs,mins,sec])'*[b_ecef_x;b_ecef_y;b_ecef_z]; % Magnetic field vector in ECI [nT]
b_b = C_bG*b_eci; % Magnetic field vector in body [nT]
T_m = cross(m_b,(b_b/(10^9))); % Magnetic torque [kN*km]

%% Aerodynamic drag with NRLMSISE00
[temp, rho] = atmosnrlmsise00(rlla(3), rlla(1), rlla(2), 2019, 151, sec,'None');
v = norm(v_b)*1000; % Velocity magnitude
v_hat = (v_b*1000)/v; % Velocity unit vector
n_x_a = dot(n_x_b_unit,v_hat); % amount of x area
n_y_a = dot(n_y_b_unit,v_hat); % amount of y area
n_z_a = dot(n_z_b_unit,v_hat); % amounf of z area
n_x_a_neg = dot(n_x_b_unit_neg,v_hat); % amount of x area
n_y_a_neg = dot(n_y_b_unit_neg,v_hat); % amount of y area
n_z_a_neg = dot(n_z_b_unit_neg,v_hat); % amounf of z area
% Checking for direction and applicable area
% Rho vectos from the center of mass to center of surface
% Cube
rho_cube_x_a = [1;0;0]*a_cube;
rho_cube_y_a = [0;1;0]*a_cube;
rho_cube_z_a = [0;0;1]*a_cube;
rho_cube_x_a_neg = [-1;0;0]*a_cube;
rho_cube_y_a_neg = [0;-1;0]*a_cube;
rho_cube_z_a_neg = [0;0;-1]*a_cube;

% Solar panel 1
rho_solar_z_a_1 = [0;2.5;0]*a_solar;
rho_solar_z_a_neg_1 = [0;2.5;0]*a_solar;

% Solar panel 2
rho_solar_z_a_2 = [0;-2.5;0]*a_solar;
rho_solar_z_a_neg_2 = [0;-2.5;0]*a_solar;

% Sensor
rho_sensor_x_a = [0.125;0;2]*a_sensor_side;
rho_sensor_y_a = [0;0.125;2]*a_sensor_side;
rho_sensor_z_a = [0;0;3]*a_sensor_bottom;
rho_sensor_x_a_neg = [-0.125;0;2]*a_sensor_side;
rho_sensor_y_a_neg = [0;-0.125;2]*a_sensor_side;

if n_x_a<0
    n_x_a = 0;
    rho_cube_x_a = [0;0;0];
    rho_sensor_x_a = [0;0;0];
end
if n_y_a<0
    n_y_a = 0;
    rho_cube_y_a = [0;0;0];
    rho_sensor_y_a = [0;0;0];
end
if n_z_a<0
    n_z_a = 0;
    rho_cube_z_a = [0;0;0];
    rho_sensor_z_a = [0;0;0];
    rho_solar_z_a_1 = [0;0;0];
    rho_solar_z_a_2 = [0;0;0];

end
if n_x_a_neg<0
    n_x_a_neg = 0;
    rho_cube_x_a_neg = [0;0;0];
    rho_sensor_x_a_neg = [0;0;0];
end
if n_y_a_neg<0
    n_y_a_neg = 0;
    rho_cube_y_a_neg = [0;0;0];
    rho_sensor_y_a_neg = [0;0;0];
end
if n_z_a_neg<0
    n_z_a_neg = 0;
    rho_cube_z_a_neg = [0;0;0];
    rho_solar_z_a_neg_1 = [0;0;0];
    rho_solar_z_a_neg_2 = [0;0;0];
end
rho_tot_a = (rho_cube_x_a*n_x_a+ ...
    rho_cube_y_a*n_y_a + ...
    rho_cube_z_a*n_z_a + ...
    rho_cube_x_a_neg*n_x_a_neg + ...
    rho_cube_y_a_neg*n_y_a_neg + ...
    rho_cube_z_a_neg*n_z_a_neg + ...
    rho_solar_z_a_1*n_z_a + ...
    rho_solar_z_a_neg_1*n_z_a_neg + ...
    rho_solar_z_a_2*n_z_a + ...
    rho_solar_z_a_neg_2*n_z_a_neg + ...
    rho_sensor_x_a*n_x_a + ...
    rho_sensor_y_a*n_y_a + ...
    rho_sensor_z_a*n_z_a + ...
    rho_sensor_x_a_neg*n_x_a_neg + ...
    rho_sensor_y_a_neg*n_y_a_neg);
a_ram = a_x*(n_x_a+n_x_a_neg)+a_y*(n_y_a+n_y_a_neg)+a_z*(n_z_a+n_z_a_neg); % Ram area [m^2]
c_ap = (rho_tot_a/a_ram); % Center of aerodynamic pressure 
F_ad = -0.5*rho(6)*a_ram*v^2*C_D*v_hat; % Aerodynamic drag force [(kg*km)/s^2]
T_ad = cross(c_ap,F_ad); % Aerodynamic drag torque [kN*km]

%% Total disturbance torque
T = T_ad + T_sp + T_gg + T_m; % Total torque [kN*km]
 % Quaternion errors
    q_bG = [eta_bG;eps_bG]; % Quaternion in MATLAB form
    q_bLVLH = [eta_bLVLH;eps_bLVLH]; % Quaternion in MATLAB form
    q_e_bLVLH = quatmultiply(quatconj(q_c_bLVLH),q_bLVLH'); % Quaternion error
    T_c = -K_p*[q_e_bLVLH(2);q_e_bLVLH(3);q_e_bLVLH(4)]-K_d*omega_bLVLH;
    deps_bG = 0.5*(q_bG(1)*eye(3)+cross_x(q_bG(2),q_bG(3),q_bG(4)))*omega_bG_sc; % Derivative of quaternion vector
    deta_bG = -0.5*[q_bG(2);q_bG(3);q_bG(4)]'*omega_bG_sc; % Derivative of quaternion scalar
    deps_bLVLH = 0.5*(q_e_bLVLH(1)*eye(3)+cross_x(q_e_bLVLH(2),q_e_bLVLH(3),q_e_bLVLH(4)))*omega_bLVLH; % Quaternion vector derivative
    deta_bLVLH = -0.5*[q_e_bLVLH(2);q_e_bLVLH(3);q_e_bLVLH(4)]'*omega_bLVLH; % Quaternion scalar derivative
    domega_bG_sc = (I_sc^-1)*(T-cross_x(omega_bG_sc(1),omega_bG_sc(2),omega_bG_sc(3))*I_sc*omega_bG_sc+T_c); % Derivative of angular velocity [rad/s^2]
    domega_bG_w = (I_w^-1)*(-cross_x(omega_bG_sc(1),omega_bG_sc(2),omega_bG_sc(3))*I_w*omega_bG_w-T_c); % Derivative of angular velocity [rad/s^2]
    dstatedt = [dx;dy;dz;ddx;ddy;ddz;deps_bG;deta_bG;deps_bLVLH;deta_bLVLH;domega_bG_sc;domega_bG_w]; % Return state vectors
end

end
