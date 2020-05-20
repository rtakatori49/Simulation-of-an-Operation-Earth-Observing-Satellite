function [M] = RBMotionSim(psi, theta, phi)

axis square, xlim([-1 1]), ylim([-1 1]), zlim([-1 1])
view(125, 50);
grid on

light('position', [1,1,0], 'HandleVisibility', 'off');
lighting gouraud

% Define initial YAW angle
psi_0 = psi(1);
% Define initial PITCH angle
theta_0 = theta(1);
% Define initial ROLL angle
phi_0 = phi(1);

% Draw Body Graphical Components
%m = 1500;
r = .5;
l = 1;

b1_vec = arrow([0 1], [0 0], [0 0], 'facecolor', 'red');
b2_vec = arrow([0 0], [0 1], [0 0], 'facecolor', 'red');
b3_vec = arrow([0 0], [0 0], [0 1], 'facecolor', 'red');
rb = cylin(r,r,l, 120, 'top', 'facecolor',[.1 .1 .1]);
translate(rb, 0,0,-l/2);

% Draw inertial graphical componenets
i1_vec = arrow([0 1], [0 0], [0 0], 'facecolor', 'blue');
i2_vec = arrow([0 0], [0 1], [0 0], 'facecolor', 'blue');
i3_vec = arrow([0 0], [0 0], [0 1], 'facecolor', 'blue');

% Get everything into initial position
rotate([b1_vec, b2_vec, b3_vec rb], [0 0 1]', psi_0*180/pi, [0 0 0]');
rotate([b1_vec, b2_vec, b3_vec rb], Cz(psi_0)'*[0 1 0]', theta_0*180/pi, [0 0 0]');
rotate([b1_vec, b2_vec, b3_vec rb], Cz(psi_0)'*Cy(theta_0)'*[1 0 0]', phi_0*180/pi, [0 0 0]');

% Simulate the motion

M(1) = getframe;
% find the deltas
dpsi = diff(psi);
dtheta = diff(theta);
dphi = diff(phi);

% loop through the results for the simulation
for i = 1:length(dtheta)
    rotate([b1_vec, b2_vec, b3_vec rb],...
        [0 0 1]', dpsi(i)*180/pi, [0 0 0]');
    rotate([b1_vec, b2_vec, b3_vec rb],...
        Cz(psi(i))'*[0 1 0]', dtheta(i)*180/pi, [0 0 0]');
    rotate([b1_vec, b2_vec, b3_vec rb],...
        Cz(psi(i))'*Cy(theta(i))'*[1 0 0]', dphi(i)*180/pi, [0 0 0]');
    M(i+1) = getframe;
    drawnow
end