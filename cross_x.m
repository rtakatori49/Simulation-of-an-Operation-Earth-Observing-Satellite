function [cross] = cross_x(x,y,z)
    cross = [0, -z, y;
        z, 0, -x;
        -y, x, 0];
end