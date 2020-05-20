function status = outputfcndetumble(t, dy, flag, mu,I_sc,dt_K_p,dt_K_d,q_c_bLVLH)

global Torques T_c

status = 0;

if strcmp(flag,'init')
    Torques = zeros(1,4);
    
elseif isempty(flag)
    Torques = [Torques; [T_c' t(end)]];

elseif strcmp(flag,'done')

end

end %end output function

