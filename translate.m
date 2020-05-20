function translate(h, x, y, z)
%
% Eric A. Mehiel
% Cal Poly, SLO
% Aerospace Engineering Department

for i = 1:length(h)
    if (strcmpi(get(h(i), 'type'), 'text'))
        p = get(h(i), 'position');
        set(h(i), 'position', p + [x y z]);
    else
        xdata = get(h(i), 'xdata') + x;
        ydata = get(h(i), 'ydata') + y;
        zdata = get(h(i), 'zdata') + z;
        set(h(i), 'xdata', xdata, 'ydata', ydata, 'zdata', zdata);
    end

end