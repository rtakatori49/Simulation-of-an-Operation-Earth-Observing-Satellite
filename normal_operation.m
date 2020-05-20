function [dstatedt] = normal_operation(t, state, mu,a_x,a_y,a_z,a_solar,a_cube,a_sensor_side,a_sensor_bottom,sp,sun,C_D, mdp,I_sc,I_w,K_p,K_d,q_c_bLVLH)
global Torques T T_ad T_sp T_gg T_m T_c
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

