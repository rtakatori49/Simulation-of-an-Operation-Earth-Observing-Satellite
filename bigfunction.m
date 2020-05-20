%% Big function
% This function is used for ode45 for homework 6 and 7. This accounts for
% the two body propagation, quaternion, body transformation, and
% disturbance torques.

function [ dstatedt ] = bigfunction( t, state, mu, I,rho_c,a_x,a_y,a_z,a_solar,a_cube,a_sensor_side,a_sensor_bottom,sp,sun,C_D, mdp)
% Gobale Variables
global Torques T T_ad T_sp T_gg T_m T_ad_c
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
epsilon = [state(7); state(8); state(9)]; % Quaternion vector
eta = state(10); % Quaternion scalar

% LVLH to body
epsilon_bLVLH = [state(14); state(15); state(16)];
eta_bLVLH = state(17);
C_bLVLH = quat2dcm([eta_bLVLH;epsilon_bLVLH]');

% ECI to LVLH
z_hat = -reci/r; % Z LVLH
y_hat = -cross(reci,veci)/norm(cross(reci,veci)); % Y LVLH
x_hat = cross(y_hat,z_hat); % X LVLH
C_LVLHG = [x_hat,y_hat,z_hat]'; % ECI to LVLH DCM 

% ECI to body matrix
C_bG = C_bLVLH*C_LVLHG; % ECI to body DCM

% Angular velocities in various frames
omega_bG = [state(11);state(12);state(13)]; % Angular velocity [rad/s]
omega_LVLHG = C_bG*(cross(reci,veci)/r^2);
omega_bLVLH = omega_bG-omega_LVLHG;

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

% Rho vectos from the center of mass to center of surface
% Cube
rho_cube_x = [1;0;0];
rho_cube_y = [0;1;0];
rho_cube_z = [0;0;1];
rho_cube_x_neg = [-1;0;0];
rho_cube_y_neg = [0;-1;0];
rho_cube_z_neg = [0;0;-1];

% Solar panel 1
rho_solar_z_1 = [0;2.5;0];
rho_solar_z_neg_1 = [0;2.5;0];

% Solar panel 2
rho_solar_z_2 = [0;-2.5;0];
rho_solar_z_neg_2 = [0;-2.5;0];

% Sensor
rho_sensor_x = [0.125;0;1.5];
rho_sensor_y = [0;0.125;1.5];
rho_sensor_z = [0;0;2];
rho_sensor_x_neg = [-0.125;0;1.5];
rho_sensor_y_neg = [0;-0.125;1.5];

rho_tot = (rho_cube_x + ...
    rho_cube_y + ...
    rho_cube_z + ...
    rho_cube_x_neg + ...
    rho_cube_y_neg + ...
    rho_cube_z_neg + ...
    rho_solar_z_1 + ...
    rho_solar_z_neg_1 + ...
    rho_solar_z_2 + ...
    rho_solar_z_neg_2 + ...
    rho_sensor_x + ...
    rho_sensor_y + ...
    rho_sensor_z + ...
    rho_sensor_x_neg + ...
    rho_sensor_y_neg)/1000;


%% Atmospheric drag Simplified
v = norm(v_b); % Velocity magnitude
v_hat = v_b/v; % Velocity unit vector
n_x_a = dot(n_x_b_unit,v_hat); % amount of x area
n_y_a = dot(n_y_b_unit,v_hat); % amount of y area
n_z_a = dot(n_z_b_unit,v_hat); % amounf of z area
n_x_a_neg = dot(n_x_b_unit_neg,v_hat); % amount of x area
n_y_a_neg = dot(n_y_b_unit_neg,v_hat); % amount of y area
n_z_a_neg = dot(n_z_b_unit_neg,v_hat); % amounf of z area
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
rho_sensor_x_a = [0.125;0;1.5]*a_sensor_side;
rho_sensor_y_a = [0;0.125;1.5]*a_sensor_side;
rho_sensor_z_a = [0;0;2]*a_sensor_bottom;
rho_sensor_x_a_neg = [-0.125;0;1.5]*a_sensor_side;
rho_sensor_y_a_neg = [0;-0.125;1.5]*a_sensor_side;
% Checking for direction and applicable area
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
c_ap = (rho_tot_a/a_ram)/1000; % Center of aerodynamic pressure 
F_ad = -0.5*rho_c*a_ram*v^2*C_D*v_hat*1000; % Aerodynamic drag force [(kg*km)/s^2]
T_ad_c = cross(c_ap,F_ad); % Aerodynamic drag torque [kN*km]

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
rho_sensor_x_s = [0.125;0;1.5]*a_sensor_side;
rho_sensor_y_s = [0;0.125;1.5]*a_sensor_side;
rho_sensor_z_s = [0;0;2]*a_sensor_bottom;
rho_sensor_x_s_neg = [-0.125;0;1.5]*a_sensor_side;
rho_sensor_y_s_neg = [0;-0.125;1.5]*a_sensor_side;

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
c_sp = (rho_tot_s/a_wet)/1000; % Center of solar pressure
F_s = -(sp*a_wet*sun_b)/1000; % Force due to sun [kN]
T_sp = cross(c_sp,F_s); % Solar pressure torque [kN*km]

%% Gravity gradient
T_gg = ((3*mu)/norm(r_b)^5)*cross_x(r_b(1),r_b(2),r_b(3))*I*r_b; % Gravity gradient torque [N*km]

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
T_m = cross(m_b,(b_b/(10^9)))/(10^6); % Magnetic torque [kN*km]

%% Aerodynamic drag with NRLMSISE00
[temp, rho] = atmosnrlmsise00(rlla(3), rlla(1), rlla(2), 2019, 151, sec,'None');
v = norm(v_b); % Velocity magnitude
v_hat = v_b/v; % Velocity unit vector
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
rho_sensor_x_a = [0.125;0;1.5]*a_sensor_side;
rho_sensor_y_a = [0;0.125;1.5]*a_sensor_side;
rho_sensor_z_a = [0;0;2]*a_sensor_bottom;
rho_sensor_x_a_neg = [-0.125;0;1.5]*a_sensor_side;
rho_sensor_y_a_neg = [0;-0.125;1.5]*a_sensor_side;

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
c_ap = (rho_tot_a/a_ram)/1000; % Center of aerodynamic pressure 
F_ad = -0.5*rho(6)*a_ram*v^2*C_D*v_hat*1000; % Aerodynamic drag force [(kg*km)/s^2]
T_ad = cross(c_ap,F_ad); % Aerodynamic drag torque [kN*km]

%% Total disturbance torque
T = T_ad + T_sp + T_gg + T_m; % Total torque [kN*km]

%% Quaternion derivatives
deps = 0.5*(eta*eye(3)+cross_x(epsilon(1),epsilon(2),epsilon(3)))*omega_bG; % Quaternion vector derivative
deta = -0.5*epsilon'*omega_bG; % Quaternion scalar derivative
deps_bLVLH = 0.5*(eta_bLVLH*eye(3)+cross_x(epsilon_bLVLH(1),epsilon_bLVLH(2),epsilon_bLVLH(3)))*omega_bLVLH; % Quaternion vector derivative
deta_bLVLH = -0.5*epsilon_bLVLH'*omega_bLVLH; % Quaternion scalar derivative
domega= (I^-1)*(T-(cross_x(omega_bG(1),omega_bG(2),omega_bG(3))*I*omega_bG)); % Angular acceleration [rad/s^2]

%% Return state
dstatedt = [dx;dy;dz;ddx;ddy;ddz;deps;deta;domega;deps_bLVLH;deta_bLVLH]; % Return state vector for next step

end

