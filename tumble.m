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

