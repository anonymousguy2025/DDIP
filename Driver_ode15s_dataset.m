% =========================================================================
% Author: YUKI
% Last Edited: 2025-03-02
% Description: This script sets up a 3-link continuum robot model, defines
% its material and geometric properties, solves for the static configuration,
% and prepares initial states for dynamics or learning-based simulation.
% =========================================================================

%% Initial conditions
startup;
clc
load('3LinkRobot.mat','S1','q0')
qd0 = zeros(S1.ndof,1);

% Number of samples to generate
num_samples = 1000;

% Set file path for saving the dataset
filepath = 'D:\yyq\DIDP_matlab2\DIDP_matlab/data_0.0012/';  % Modify this to your desired folder path
% Create the folder if it does not exist
if ~exist(filepath, 'dir')
    mkdir(filepath);
end


%% Generate random theta samples (2 joints, 4 time steps each)
theta_samples = zeros(num_samples, 2, 4);

for i = 1:num_samples
    if mod(i, 500) == 0  % Check if i is divisible by 500
        rng('shuffle'); % Shuffle the random seed based on the current time
    end
    for j = 1:4
        theta_samples(i, 1, j) = -pi + (2 * pi) * rand();               % joint 1: [-π, π]
        theta_samples(i, 2, j) = -pi/2 + (pi/4 - (-pi/2)) * rand();     % joint 2: [-π/2, π/4]
    end
end



%% Main simulation loop
k = 1;

for i = 1:num_samples
    T = 2; % Duration for actuation
    theta = squeeze(theta_samples(i, :, :)); % Extract theta (2x4)
    uqt = @(t) ActuationInput(t,T,theta); % Create actuation input function

    % Joint angle bounds
    lower_bound1 = -pi;
    upper_bound1 = pi;
    lower_bound2 = -pi/4;
    upper_bound2 = pi/2;


    skip_save = false;
    t_check = linspace(0, T, 1500); % Time steps to verify joint limits

    % Validate theta ranges over the full duration
    for tt = 1:length(t_check)
        uqt_val = uqt(t_check(tt));
        theta1 = uqt_val(1,1);
        theta2 = uqt_val(2,1);

        if theta1 < lower_bound1 || theta1 > upper_bound1
            fprintf('Sample %d at t=%.3f: theta1 is out of range (%.3f), skipping save.\n', i, t_check(tt), theta1);
            skip_save = true;
            break;
        end

        if theta2 < lower_bound2 || theta2 > upper_bound2
            fprintf('Sample %d at t=%.3f: theta2 is out of range (%.3f), skipping save.\n', i, t_check(tt), theta2);
            skip_save = true;
            break;
        end
    end

    % If the sample is out of range, skip it and do not execute the following code
    if skip_save
        fprintf('Sample %d does not meet the criteria, skipping.\n', i);
        continue; % **Skip the current sample and proceed to the next i**
    end

    fprintf('Sample %d data is within range and can be saved.\n', i);

    % Compute initial state
    qqdqddj0 = uqt(0);
    q0(S1.i_jactq)  = qqdqddj0(:,1); %replacing q_k and qdot_k
    qd0(S1.i_jactq) = qqdqddj0(:,2);
    qqd0 = [q0;qd0];

    qddu0 = FDU_IDM(S1,0,q0,qd0,uqt);
    qdd0 = qddu0(1:S1.ndof);
    u0 = qddu0(S1.ndof+1:end);
    qqdu0 = [q0;qd0;u0];
    qdqddud0 = [qd0;qdd0;u0*0]; %ud0 never used

    % ODE parameters
    tmax = T;
    dt = 0.001;

    %% method 3.3 with Jacobian, multiple function call
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'Jacobian',@(t,qqd) ODEJacobian(S1,t,qqd,FD_SoRoSim(S1,t,qqd,uqt),uqt));% may not be the most efficient implementation, qdd may be computed twice and derivatives require additional computations
    tic
    % Solve system dynamics
    [t1,qqd1] = ode15s(@(t,qqd) derivatives_SoRoSim(S1,t,qqd,uqt),0:dt:tmax,qqd0,options);
    toc

    t=t1;
    qqd=qqd1;
    q = qqd(:,1:20);
    qd = qqd(:,21:40);

    num_samples_q = size(q, 1);
    all_positions = zeros(num_samples_q, 21);
    all_velocities = zeros(num_samples_q, 21);

    %% positions, velocities
    for ii = 1:num_samples_q

        current_q = q(ii, :);
        current_qd = qd(ii, :);

        [positions, velocities] = PositionVelocity(S1, current_q', current_qd', 0); % 转置为列向量

        all_positions(ii, :) = positions(1:21)'; %positions 的前 3 个分量是 x, y, z
        all_velocities(ii, :) = velocities(1:21)'; % velocities 的前 3 个分量是 vx, vy, vz
    end

    %% saving to HDF5
    h5_filename = [filepath,'DynamicsSolution_SoRoSim_', num2str(k), '.h5'];

    % 统一 ChunkSize 行数为 100，列数根据数据自动匹配
    fixed_chunk_size = @(data) [min(100, size(data, 1)), size(data, 2)];

    % 创建 HDF5 并写入数据（t, qqd, theta, T）
    h5create(h5_filename, '/t', size(t), 'Datatype', 'double', 'ChunkSize', fixed_chunk_size(t), 'Deflate', 9);
    h5write(h5_filename, '/t', t);

    h5create(h5_filename, '/qqd', size(qqd), 'Datatype', 'double', 'ChunkSize', fixed_chunk_size(qqd), 'Deflate', 9);
    h5write(h5_filename, '/qqd', qqd);

    h5create(h5_filename, '/theta', size(theta), 'Datatype', 'double', 'ChunkSize', fixed_chunk_size(theta), 'Deflate', 9);
    h5write(h5_filename, '/theta', theta);

    h5create(h5_filename, '/T', size(T), 'Datatype', 'double', 'ChunkSize', fixed_chunk_size(T), 'Deflate', 9);
    h5write(h5_filename, '/T', T);

    % 使用固定的 ChunkSize = [100, 列数]
    chunk_size_positions = fixed_chunk_size(all_positions);
    chunk_size_velocities = fixed_chunk_size(all_velocities);

    h5create(h5_filename, '/all_positions', size(all_positions), 'Datatype', 'double', 'ChunkSize', chunk_size_positions, 'Deflate', 9);
    h5write(h5_filename, '/all_positions', all_positions);

    h5create(h5_filename, '/all_velocities', size(all_velocities), 'Datatype', 'double', 'ChunkSize', chunk_size_velocities, 'Deflate', 9);
    h5write(h5_filename, '/all_velocities', all_velocities);
    out_folder = fullfile(filepath, 'videos');
    if ~exist(out_folder, 'dir')
        mkdir(out_folder);
    end
    if k < 11 || mod(k, 500) == 0
        plotqqd_train(S1,t,qqd,k,out_folder)
    end
    k = k + 1;
end