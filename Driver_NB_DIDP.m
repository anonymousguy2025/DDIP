% =========================================================================
% Author: YUKI
% Last Edited: 2025-04-15
% Description: Evaluate prediction accuracy of action trajectories via SoRoSim.
%              Computes Euclidean distance to ground truth goal, tracks min.
% =========================================================================

%% Initial conditions
startup
clc
load('3LinkRobot','S1','q0')

%% ========================= File and variable initialization =========================
folderPath = 'D:\yyq\ThreeLinkHybrid\500\h5_file\';% Modify this to your desired folder path
files = dir(fullfile(folderPath, 'predict_*.h5'));
total_count = numel(files) - 1;

% Accuracy counters
correct_count_5cm = 0;
correct_count_10cm = 0;
correct_count_2cm = 0;
correct_count_1cm = 0;
ave_distance = 0;

%% ========================= Loop Over Each Prediction File =========================
for i = 0:total_count
    % -------- Load prediction file --------
    fileName = sprintf('predict_%d.h5', i);
    filePath = fullfile(folderPath, fileName);
    fprintf('Processing: %s\n', filePath);
    h5disp(filePath);

    % Read dataset
    action_gt = h5read(filePath, '/action_gt');
    action_pred = h5read(filePath, '/action_pred');
    target_pos = h5read(filePath, '/target_pos');

    % Use only first two rows (2 joints)
    squeeze(action_gt(1, :, :))
    action_gt = action_gt(1:2,:);
    squeeze(action_pred(1, :, :))
    action_pred = action_pred(1:2,:);


    % -------- Trajectory and smoothing setup --------
    T_NB = 0.5;% Control duration
    dt = 0.001;% Time step
    nt = T_NB / dt+1; % nt should be > 1

    xyz_gt_pos = target_pos;
    xyz_gt_pos = xyz_gt_pos(:); % Convert to column vector

    % 目标时间向量（新的时间步长为 0.0005）
    dt = 0.001;
    t_new = 0:dt:T_NB;
    t_old = 0:dt:T_NB-dt;

    T_smooth = 0.05;
    ts = 0:dt:T_smooth;

    % Cosine window for smooth startup
    alpha = ones(1, length(t_new));
    alpha(1:length(ts)) = 0.5 * (1 - cos(pi * ts / T_smooth));

    action_pred_interp = zeros(size(action_pred, 1), length(t_new));
    for k = 1:size(action_pred, 1)
        action_pred_interp(k, :) = interp1(t_old, action_pred(k, :), t_new, 'linear');
    end

    action_pred = action_pred_interp;

    % Smooth and apply window
    action_pred(1,:) = smooth(action_pred(1,:),20).*alpha';
    action_pred(2,:) = smooth(action_pred(2,:),20).*alpha';

    theta_NB = action_pred';

    % ===================== Initial state setup ======================
    qd0 = zeros(S1.ndof,1);
    qqd0 = [q0;qd0];
    thetadd0 = [0 0];
    qddu0 = FDU_SoRoSim(S1,0,qqd0,thetadd0);
    qdd0 = qddu0(1:S1.ndof);
    u0 = qddu0(S1.ndof+1:S1.ndof+S1.nact);

    q_pre = q0;
    u_pre = u0;
    qu_pre = [q_pre;u_pre];
    qd_pre = qd0;
    qdd_pre= qdd0;
    qqdqdd_pre = [q_pre;qd_pre;qdd_pre];

    % Newmark-beta integration parameters
    beta = 0.25;
    gamma = 0.5;

    tmax = T_NB;
    t=0:dt:tmax;
    if t(end)~=tmax
        t=[t,tmax];
    end
    nt = length(t);
    qqdqddu = zeros(nt, 3*S1.ndof+S1.nact);

    % Solver setup
    options = optimoptions('fsolve','Jacobian','on','Display','final','Algorithm','trust-region-dogleg','MaxIterations',1000); % Set Display to 'final'/iter to see if solved. 'off' for none.
    qu_next_guess = qu_pre;
    exflg = zeros(nt,1);
    qqdqddu(1,:) = [qqdqdd_pre' u0'];

    tic

    % Track minimum distance to target
    min_distance = inf;
    best_position = [];

    % Failure detection flags
    first_fail_j = -1;
    consecutive_fail_count = 0;
    is_failed = false;

    % ===================== Simulation Loop ======================
    for j=2:nt
        [qu,~,exflg(j)] = fsolve(@(qu) IDResidueJacobian(S1,t(j),qu,dt,qqdqdd_pre,theta_NB(j-1,:)), qu_next_guess, options);

        % Failure handling
        if exflg(j) <= 0
            consecutive_fail_count = consecutive_fail_count + 1;
            if first_fail_j == -1
                first_fail_j = j;
            end
        end

        % Check if there have been 10 consecutive failures since the first failure
        if first_fail_j > 0 && (j - first_fail_j + 1) >= 10
            fprintf('⚠️  File %d failed to solve 10 consecutive steps starting from %.3f (index j = %d), breaking loop.\n', ...
                i, t(first_fail_j), first_fail_j);

            % Write to fail log
            folderPath = regexprep(folderPath, '[\\/]+$', '');
            [parentFolder1, ~, ~] = fileparts(folderPath);
            [~, log_prefix, ~] = fileparts(parentFolder1);
            log_file = ['fail_log_' log_prefix '.txt'];
            fileID = fopen(log_file, 'a');
            fprintf(fileID, '[%s] File %d failed 10 consecutive times from %.4f s (index %d)\n', ...
                datestr(now, 'yyyy-mm-dd HH:MM:SS'), i, t(first_fail_j), first_fail_j);
            fclose(fileID);
            is_failed = true;  % Mark the file as failed
            break;
        end

        % If successful step, update state using Newmark-beta
        if exflg(j) > 0
            q = qu(1:S1.ndof);
            u = qu(S1.ndof+1:S1.ndof+S1.nact);
            qd  = gamma/(beta*dt)*(q-q_pre)+(1-gamma/beta)*qd_pre+dt*(1-gamma/(2*beta))*qdd_pre;
            qdd = 1/(beta*dt^2)*(q-q_pre)-1/(beta*dt)*qd_pre+(1-1/(2*beta))*qdd_pre;

            qqdqddu(j,:) = [q;qd;qdd;u]';

            % Update previous states
            q_pre = q;
            u_pre = u;

            qu_pre = [q_pre;u_pre];
            qd_pre = qd;
            qdd_pre = qdd;
            qqdqdd_pre = [q_pre;qd_pre;qdd_pre];
            qu_next_guess = qu_pre;

            % -------- Update closest prediction --------
            [positions, velocities] = PositionVelocity(S1, q', qd', 0);
            xyz_pred_pos = positions(19:21);
            xyz_pred_pos = xyz_pred_pos(:);

            distance_tmp = norm(xyz_gt_pos - xyz_pred_pos, 2);
            if distance_tmp < min_distance
                min_distance = distance_tmp;
                best_position = xyz_pred_pos;
            end
            t(j) % Debug display
        end
    end
    toc

    % ===================== Accuracy Evaluation ======================
    distance = min_distance;
    ave_distance = ave_distance + distance
    fprintf('Euclidean distance between xyz_gt_pos and xyz_pred_pos: %.6f\n', distance);
    if distance <= 0.1
        correct_count_10cm = correct_count_10cm + 1;
    end

    if distance <= 0.05
        correct_count_5cm = correct_count_5cm + 1;
    end

    if distance <= 0.02
        correct_count_2cm = correct_count_2cm + 1;
    end

    if distance <= 0.01
        correct_count_1cm = correct_count_1cm + 1;
    end

    % ===================== Final Statistics ======================
    if i == total_count
        accuracy_10cm = (correct_count_10cm / (total_count+1)) * 100;
        accuracy_5cm = (correct_count_5cm / (total_count+1)) * 100;
        accuracy_2cm = (correct_count_2cm / (total_count+1)) * 100;
        accuracy_1cm = (correct_count_1cm / (total_count+1)) * 100;
        ave_distance = ave_distance / (total_count+1);

        % Append results to output file
        fileID = fopen('results.txt', 'a');
        fprintf(fileID, '\n==== Result on 2dof %s ====\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));

        fprintf(fileID, 'Final 10cm Accuracy: %.2f%%\n', accuracy_10cm);
        fprintf(fileID, 'Final 5cm Accuracy: %.2f%%\n', accuracy_5cm);
        fprintf(fileID, 'Final 2cm Accuracy: %.2f%%\n', accuracy_2cm);
        fprintf(fileID, 'Final 1cm Accuracy: %.2f%%\n', accuracy_1cm);
        fprintf(fileID, 'ave_distance: %.7f\n', ave_distance);

        fclose(fileID);

        fprintf('Final 10cm Accuracy: %.2f%%\n', accuracy_10cm);
        fprintf('Final 5cm Accuracy: %.2f%%\n', accuracy_5cm);
        fprintf('Final 2cm Accuracy: %.2f%%\n', accuracy_2cm);
        fprintf('Final 1cm Accuracy: %.2f%%\n', accuracy_1cm);
        fprintf('ave_distance: %.7f\n', ave_distance);
    end
    %plot(t,exflg)
    % ===================== Save result figure ======================
    if ~is_failed
        out_folder = fullfile(fileparts(folderPath), 'pred_videos');
        if ~exist(out_folder, 'dir'), mkdir(out_folder); end

        plotqqd_test(S1, t, qqdqddu, i, xyz_gt_pos, xyz_pred_pos, out_folder);
    else
        fprintf('⚠️ Skip saving video for file %d due to simulation failure.\n', i);
    end

end
