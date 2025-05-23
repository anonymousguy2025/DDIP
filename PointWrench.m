function F = PointWrench(t)
F = [0.1*cos(t)^2,-0.2,0.3*sin(t),-2+t/2.5,4-t,-2]';
if t>=5
    F(2) = 0.2;
    F(5) = -6+t;
    F(6) = 2;
end
% F = [0,0,0,0,0,0]';
end