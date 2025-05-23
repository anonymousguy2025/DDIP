function q_t = HermiteInterpolate(tt, t, q, qd)
% Find the interval in which t falls
idx = find(tt >= t, 1, 'last');

if idx == length(t)
    idx = length(t) - 1;
end

% Extract the relevant times and values
t1 = t(idx);
t2 = t(idx+1);
q1 = q(:,idx);
q2 = q(:,idx+1);
qd1 = qd(:,idx);
qd2 = qd(:,idx+1);

% Number of elements in the vector
n = length(q1);

% Initialize the result vector
q_t = zeros(n, 1);

% Normalize the time
duration = t2 - t1;
tau = (tt - t1) / duration;  % Normalized time

% Construct the matrix for the system of equations
M = [1, 0, 0, 0;
     0, 1, 0, 0;
     1, 1, 1, 1;
     0, 1, 2, 3];

% Loop over each element in the vector
for i = 1:n
    % Construct the vector for the known values
    b = [q1(i); qd1(i); q2(i); qd2(i)];
    
    % Solve the system for the coefficients
    a = M \ b;
    
    % Evaluate the polynomial at normalized time tau
    q_t(i) = a(1) + a(2)*tau + a(3)*tau^2 + a(4)*tau^3;
end

end