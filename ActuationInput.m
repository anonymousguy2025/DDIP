% =========================================================================
% Author: YUKI
% Last Edited: 2025-02-11
% Description: This function computes the joint positions, velocities, and
% accelerations of a robot's actuated joints over time using polynomial
% interpolation for smooth motion. The function takes in the current time
% (t), total time (T), and a matrix of joint motion coefficients (theta),
% and outputs the joint configurations and their derivatives.
% =========================================================================

function qjt = ActuationInput(t,T,theta)
%theta is parameterized value of theta1 and theta2

% Time Intervals & Polynomial Constants
Ti = T/4;
c1 = (Ti)^3*(Ti-T/2)*(Ti-3 *T/4)*(Ti-T);
Ti = T/2;
c2= (Ti)^3*(Ti-T/4)*(Ti-3*T/4)*(Ti-T);
Ti = 3*T/4;
c3= (Ti)^3*(Ti-T/4)*(Ti-T/2)*(Ti-T);
Ti = T;
c4= (Ti)^3*(Ti-T/4)*(Ti-T/2)*(Ti-3*T/4);

% Polynomial Functions for Position, Velocity, and Acceleration
f1 = t^3*(t-T/2)*(t-3*T/4)*(t-T)/c1;
f2 = t^3*(t-T/4)*(t-3*T/4)*(t-T)/c2;
f3 = t^3*(t-T/4)*(t-T/2)*(t-T)/c3;
f4 = t^3*(t-T/4)*(t-T/2)*(t-3*T/4)/c4;

% First derivatives for velocity (polynomial velocity equations)
fd1 = -(t^2*(9*T^3 - 52*T^2*t + 90*T*t^2 - 48*t^3))/(8*c1);
fd2 = -(t^2*(9*T^3 - 76*T^2*t + 160*T*t^2 - 96*t^3))/(16*c2);
fd3 = -(t^2*(3*T^3 - 28*T^2*t + 70*T*t^2 - 48*t^3))/(8*c3);
fd4 = -(t^2*(9*T^3 - 88*T^2*t + 240*T*t^2 - 192*t^3))/(32*c4);

% Second derivatives for acceleration (polynomial acceleration equations)
fdd1 = -(3*t*(3*T^3 - 26*T^2*t + 60*T*t^2 - 40*t^3))/(4*c1);
fdd2 = -(t*(9*T^3 - 114*T^2*t + 320*T*t^2 - 240*t^3))/(8*c2);
fdd3 = -(t*(3*T^3 - 42*T^2*t + 140*T*t^2 - 120*t^3))/(4*c3);
fdd4 = -(3*t*(3*T^3 - 44*T^2*t + 160*T*t^2 - 160*t^3))/(16*c4);

% Joint Configuration Calculation
qjt = zeros(2,3);
qjt(:,1) = [f1*theta(1,1)+f2*theta(1,2)+f3*theta(1,3)+f4*theta(1,4);...
    f1*theta(2,1)+f2*theta(2,2)+f3*theta(2,3)+f4*theta(2,4)];
qjt(:,2) = [fd1*theta(1,1)+fd2*theta(1,2)+fd3*theta(1,3)+fd4*theta(1,4);...
    fd1*theta(2,1)+fd2*theta(2,2)+fd3*theta(2,3)+fd4*theta(2,4)];
qjt(:,3) = [fdd1*theta(1,1)+fdd2*theta(1,2)+fdd3*theta(1,3)+fdd4*theta(1,4);...
    fdd1*theta(2,1)+fdd2*theta(2,2)+fdd3*theta(2,3)+fdd4*theta(2,4)];

% Final Check: Ensure the motion stops if the time exceeds the total time
% the final joint positions with zero velocity and acceleration.
if t>T
    qjt(1,1) = theta(1,4); qjt(1,2) = 0; qjt(1,3) = 0;
    qjt(2,1) = theta(2,4); qjt(2,2) = 0; qjt(2,3) = 0;
end

