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

