%% COE to RV
% Ryo Takatori

function [ rperi, vperi, reci, veci ] = coe2rv( mu, e, h, inc, OMEGA, omega, theta)
% State vectors in perifocal
rperi = ((h^2)/mu)*(1/(1+e*cosd(theta)))*[cosd(theta) sind(theta) 0];
vperi = (mu/h)*[-sind(theta) e+cosd(theta) 0];

% Perifocal to ECI
[ reci, veci ] = peri2eci( rperi, vperi, omega, OMEGA, inc );

end

