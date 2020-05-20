function status = outputfcnFinal(t, dy, flag, mu,a_x,a_y,a_z,a_solar,a_cube,a_sensor_side,a_sensor_bottom,sp,sun_ECI,C_D, mdp,I_sc_tot,I_w,K_p,K_d,q_c_bLVLH)

global Torques T T_ad T_sp T_gg T_m T_c

status = 0;

if strcmp(flag,'init')
    Torques = zeros(1,19);
    
elseif isempty(flag)
    Torques = [Torques; [T' T_ad' T_sp' T_gg' T_m' T_c' t(end)]];

elseif strcmp(flag,'done')

end

end %end output function

