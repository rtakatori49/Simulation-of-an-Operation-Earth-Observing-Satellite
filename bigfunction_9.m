function [dstatedt] = bigfunction_9(t,state,mu,I_w,K_p,K_d,q_c_bLVLH)
dx = state(4); % dx [km/s]
dy = state(5); % dy [km/s]
dz = state(6); % dz [km/s]

reci = [state(1);state(2);state(3)];
veci = [state(4);state(5);state(6)];
r = norm(reci); % Magnitude of position [km]
v = norm(veci);
% Equation of motion
ddx = -mu*state(1)/r^3; % ddx [km/s^2]
ddy = -mu*state(2)/r^3; % ddy [km/s^2]
ddz = -mu*state(3)/r^3; % ddz [km/s^2]
eps_bG = [state(7);state(8);state(9)]; % Quaternion vector

%% Quaternion and Rotation Matrix
% LVLH to body
eta_bG = state(10); % Quaternion scalar
eps_bLVLH = [state(14);state(15);state(16)];
eta_bLVLH = state(17);
C_bLVLH = quat2dcm([eta_bLVLH;eps_bLVLH]');
% ECI to LVLH
z_hat = -reci/r; % Z LVLH
y_hat = -cross(reci,veci)/norm(cross(reci,veci)); % Y LVLH
x_hat = cross(y_hat,z_hat); % X LVLH
C_LVLHG = [x_hat,y_hat,z_hat]'; % ECI to LVLH DCM 
% ECI to body
C_bG = C_bLVLH*C_LVLHG; % ECI to body DCM
% Angular velovity
    omega_bG_sc = [state(11);state(12);state(13)]; % Angular velocity [rad/s]
    omega_LVLHG = C_bG*(cross(reci,veci)/r^2);
    omega_bLVLH = omega_bG_sc-omega_LVLHG;
    omega_bG_w = [state(18);state(19);state(20)];
    % Quaternion errors
    q_bG = [eta_bG;eps_bG]; % Quaternion in MATLAB form
    q_bLVLH = [eta_bLVLH;eps_bLVLH]; % Quaternion in MATLAB form
    q_e_bLVLH = quatmultiply(quatconj(q_c_bLVLH),q_bLVLH'); % Quaternion error
    
    deps_bG = 0.5*(q_bG(1)*eye(3)+cross_x(q_bG(2),q_bG(3),q_bG(4)))*omega_bG_sc; % Derivative of quaternion vector
    deta_bG = -0.5*[q_bG(2);q_bG(3);q_bG(4)]'*omega_bG_sc; % Derivative of quaternion scalar
    deps_bLVLH = 0.5*(q_e_bLVLH(1)*eye(3)+cross_x(q_e_bLVLH(2),q_e_bLVLH(3),q_e_bLVLH(4)))*omega_bLVLH; % Quaternion vector derivative
    deta_bLVLH = -0.5*[q_e_bLVLH(2);q_e_bLVLH(3);q_e_bLVLH(4)]'*omega_bLVLH; % Quaternion scalar derivative
    domega_bG_sc = (I_sc^-1)*(-cross_x(omega_bG_sc(1),omega_bG_sc(2),omega_bG_sc(3))*I_sc*omega_bG_sc-K_p*[q_e_bG(2);q_e_bG(3);q_e_bG(4)]-K_d*omega_bLVLH); % Derivative of angular velocity [rad/s^2]
    domega_bG_w = (I_w^-1)*(-cross_x(omega_bG_sc(1),omega_bG_sc(2),omega_bG_sc(3))*I_w*omega_bG_w+K_p*[q_e_bG(2);q_e_bG(3);q_e_bG(4)]+K_d*omega_bLVLH); % Derivative of angular velocity [rad/s^2]
    dstatedt = [dx;dy;dz;ddx;ddy;ddz;deps_bG;deta_bG;domega_bG_sc;deps_bLVLH;deta_bLVLH;domega_bG_w]; % Return state vectors
end

