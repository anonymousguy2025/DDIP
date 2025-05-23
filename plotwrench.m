% Define time vector
t = linspace(0, 10, 500); % 100 points from 0 to 10

% Initialize vectors for moments and forces
Mx = zeros(size(t));
My = zeros(size(t));
Mz = zeros(size(t));
Fx = zeros(size(t));
Fy = zeros(size(t));
Fz = zeros(size(t));

% Loop through each time point to get the values from PointWrench
for i = 1:length(t)
    F = PointWrench(t(i)); % Get the 6x1 vector at time t(i)
    Mx(i) = F(1); % Extract moments
    My(i) = F(2);
    Mz(i) = F(3);
    Fx(i) = F(4); % Extract forces
    Fy(i) = F(5);
    Fz(i) = F(6);
end

% Plot Moments (Mx, My, Mz)
figure;
subplot(2, 1, 1); % First subplot for moments
plot(t, Mx, 'r', 'DisplayName', 'Mx');
hold on;
plot(t, My, 'g', 'DisplayName', 'My');
plot(t, Mz, 'b', 'DisplayName', 'Mz');
hold off;
title('Moments');
xlabel('Time (s)');
ylabel('Moment (Nm)');
legend('show');
grid on;

% Plot Forces (Fx, Fy, Fz)
subplot(2, 1, 2); % Second subplot for forces
plot(t, Fx, 'r', 'DisplayName', 'Fx');
hold on;
plot(t, Fy, 'g', 'DisplayName', 'Fy');
plot(t, Fz, 'b', 'DisplayName', 'Fz');
hold off;
title('Forces');
xlabel('Time (s)');
ylabel('Force (N)');
legend('show');
grid on;
