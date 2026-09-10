% File: MP1_solution.m
% Course: MTR 5100-Fundamentals of Robotics 
% Purpose: check workspace reachability and drive feasibility for a planar 2R robot
% Inputs: S, L1, L2, m1, m2, mp, joint/velocity limits, G, eff_g, motor ratings
% Outputs: workspace plot, reachability table, drive-feasibility table
% Author: Ali Khuzam
% Date: 10/09/2026
% Student parameter S = 5

clear; clc; close all;
%==========================================================================
%% files path
script_dir = fileparts(mfilename('fullpath'));
results_dir = fullfile(script_dir, '..', 'results');
tables_dir = fullfile(results_dir, 'tables');
figures_dir = fullfile(results_dir, 'figures');
if ~exist(tables_dir, 'dir'); mkdir(tables_dir); end
if ~exist(figures_dir, 'dir'); mkdir(figures_dir); end
%==========================================================================

%% Parameters
S = 5;
seed_val = S;

robot.L1 = 0.45; robot.L2 = 0.35; %length in (m)
robot.m1 = 2.0; robot.m2 = 1.5;   %mass in (kg)
robot.lc1 = robot.L1/2; robot.lc2 = robot.L2/2; %length in (m)
robot.I1 = (robot.m1*(robot.L1.^2))/12; robot.I2 = (robot.m2*(robot.L2.^2))/12; %inertia in (kg.m^2)
robot.mp = 0.5; %payload in (kg)

robot.q1_limit_deg = [-160, 160];   % deg
robot.q2_limit_deg = [-150, 150];   % deg
robot.q1_limit = deg2rad(robot.q1_limit_deg); % rad
robot.q2_limit = deg2rad(robot.q2_limit_deg); % rad
robot.w1_limit = 1.5; robot.w2_limit = 2; %Velocity in (rad/s)
robot.a1_limit = 2.0; robot.a2_limit = 3; %Acceleration in (rad/s^2)
robot.G = 50; %gearbox ratio
robot.eff_G = 0.85; % gearbox efficiency
robot.t_m_rate = 0.8; %rated torque in (N.m)   
robot.w_m_rate_rpm = 3000; % rated speed in (rpm)
robot.w_m_rate = robot.w_m_rate_rpm*2*pi/60; % rad/s
robot.g = 9.81; % gravity m/s^2
robot.ns = 1.5; % safety factor

% Student-specific task points (S = 5)
P1 = [0.52 + 0.004*S; 0.10 + 0.002*S];   % [x; y] in m
P2 = [0.34 + 0.002*S; 0.43 - 0.001*S];   % [x; y] in m
P3 = [0.82 + 0.003*S; 0];                % [x; y] in m
task_points = [P1, P2, P3]';             
point_labels = {'P1','P2','P3'};

%============================================
%% Joint-limited workspace sampling and FK evaluation
N = 20000;
rng(seed_val, 'twister');  
q1_samples = robot.q1_limit(1) + (robot.q1_limit(2)-robot.q1_limit(1))*rand(N,1);

q2_samples = robot.q2_limit(1) + (robot.q2_limit(2)-robot.q2_limit(1))*rand(N,1);
[x_samples, y_samples] = FK(q1_samples, q2_samples, robot);

%============================================
%% Plot workspace with task points marked
fig1 = figure;
scatter(x_samples, y_samples, 4, [0.6 0.6 0.6], 'filled'); hold on;
plot(task_points(:,1), task_points(:,2), 'r^', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
for i = 1:3
    text(task_points(i,1)+0.05, task_points(i,2), point_labels{i}, 'FontSize', 12, 'FontWeight','bold');
end
xlabel('x (m)'); ylabel('y (m)');
title('Planar 2R Manipulator: Joint-Limited Workspace (Monte Carlo Sampling)');
legend('Monte Carlo workspace samples','Task-specific target points','Location','bestoutside');
axis equal; 
axis([-1 1 -1 1]);
xticks(-1:0.2:1)
yticks(-1:0.2:1)
grid on; hold off;
% saveas(fig1, fullfile(figures_dir, 'MP1_workspace.png'));
exportgraphics(fig1, fullfile(figures_dir, 'MP1_workspace.png'), 'Resolution', 300)
%==========================================================================
%% Reachability check for each task point
r_inner = abs(robot.L1 - robot.L2);
r_outer = robot.L1 + robot.L2;
 
radial_ok      = false(3,1);
ik_geo_ok      = false(3,1);
ik_limit_ok    = false(3,1);
final_reach    = false(3,1);
q1_found_deg   = nan(3,1);
q2_found_deg   = nan(3,1);
nearest_dist   = nan(3,1);
explanation    = strings(3,1);
 
for i = 1:3
    x = task_points(i,1);
    y = task_points(i,2);
    r = sqrt(x^2 + y^2);
 
    radial_ok(i) = (r >= r_inner) && (r <= r_outer);
 
    [q1_i, q2_i, reach_geo] = IK(x, y, robot);
    ik_geo_ok(i) = reach_geo;
 
    if reach_geo
        ik_limit_ok(i) = withinJointLimits(q1_i, q2_i, robot);
        q1_found_deg(i) = rad2deg(q1_i);
        q2_found_deg(i) = rad2deg(q2_i);
    else
        ik_limit_ok(i) = false;
    end
 
    final_reach(i) = ik_geo_ok(i) && ik_limit_ok(i);
 
    d = sqrt((x_samples - x).^2 + (y_samples - y).^2);
    nearest_dist(i) = min(d);
 
    if ~radial_ok(i)
        explanation(i) = "Outside geometric annulus |L1-L2| <= r <= L1+L2 -> unreachable";
    elseif ~ik_geo_ok(i)
        explanation(i) = "Inside annulus but IK cosine condition failed -> unreachable";
    elseif ~ik_limit_ok(i)
        explanation(i) = "IK solution exists but violates joint limits -> unreachable";
    else
        explanation(i) = "Radial bound and IK/joint-limit checks both pass -> reachable";
    end
end
 
reach_table = table(point_labels', task_points(:,1), task_points(:,2), r_inner*ones(3,1), ...
    r_outer*ones(3,1), radial_ok, final_reach, q1_found_deg, q2_found_deg, ...
    nearest_dist, explanation, ...
    'VariableNames', {'Point','x_m','y_m','r_inner_m','r_outer_m','RadialBoundOK', ...
    'IK_JointLimitOK','q1_deg','q2_deg','NearestSampledDist_m','Explanation'});
 
disp(reach_table);
writetable(reach_table, fullfile(tables_dir, 'MP1_reachability_table.csv'));
%==========================================================================
%% Verification check: FK(IK(point)) round-trip for reachable points
fprintf('\n--- Kinematic Consistency Check: FK(IK(x,y)) Residual Analysis ---\n');
for i = 1:3
    if final_reach(i)
        [xv, yv] = FK(deg2rad(q1_found_deg(i)), deg2rad(q2_found_deg(i)), robot);
        resid = sqrt((xv - task_points(i,1))^2 + (yv - task_points(i,2))^2);
        fprintf('%s: positional residual = %.3e m (expected ~0 for a consistent IK solution)\n', point_labels{i}, resid);
    else
        fprintf('%s: outside the reachable workspace - no valid IK solution to verify\n', point_labels{i});
    end
end
%==========================================================================
%% Gravity torques and motor-side drive feasibility
[tau1_g, tau2_g] = gravityTorques(robot);
 
[pass1, d1] = driveFeasibility(tau1_g, robot.w1_limit, robot);
[pass2, d2] = driveFeasibility(tau2_g, robot.w2_limit, robot);
 
joint_names   = {'Joint 1','Joint 2'};
tau_g_max     = [tau1_g; tau2_g];
tau_j_req     = [d1.tau_j_req; d2.tau_j_req];
tau_m_req     = [d1.tau_m_req; d2.tau_m_req];
w_m_req_rads  = [d1.w_m_req; d2.w_m_req];
w_m_req_rpm   = w_m_req_rads*60/(2*pi);
torque_status = [d1.torque_ok; d2.torque_ok];
speed_status  = [d1.speed_ok; d2.speed_ok];
pass_status   = [pass1; pass2];
 
drive_table = table(joint_names', tau_g_max, tau_j_req, tau_m_req, w_m_req_rpm, ...
    torque_status, speed_status, pass_status, ...
    'VariableNames', {'Joint','GravityTorque_Nm','JointTorqueReq_Nm','MotorTorqueReq_Nm', ...
    'MotorSpeedReq_rpm','TorqueOK','SpeedOK','Pass'});
 
disp(drive_table);
writetable(drive_table, fullfile(tables_dir, 'MP1_drive_feasibility_table.csv'));
%==========================================================================
%% Limiting subsystem identification and end-effector selection
torque_util = tau_m_req ./ robot.t_m_rate;      % fraction of rated torque used
speed_util  = w_m_req_rads ./ robot.w_m_rate;   % fraction of rated speed used
 
util_labels = {'Joint 1 torque','Joint 1 speed','Joint 2 torque','Joint 2 speed'};
util_values = [torque_util(1), speed_util(1), torque_util(2), speed_util(2)];
[max_util, idx] = max(util_values);
 
fprintf('\n--- Engineering Summary (Student Parameter S = %d) ---\n', S);
for i = 1:3
    fprintf('%s: %s\n', point_labels{i}, explanation(i));
end
if all(pass_status)
    fprintf('Drivetrain (motor and gearbox) satisfies torque and speed requirements for both joints.\n');
else
    fprintf('Drivetrain does not satisfy requirements for at least one joint; refer to the drive-feasibility table for details.\n');
end
fprintf('Limiting subsystem: %s, operating at %.1f%% of rated capacity (smallest margin among all checks performed).\n', ...
    util_labels{idx}, max_util*100);
fprintf(['Recommended end effector: two-finger parallel gripper, suitable for a %.2f kg payload ' ...
    'in a low-complexity pick-and-place application.\n'], robot.mp);

%============================================
%% FUNCTIONS
function [x, y] = FK(q1, q2, robot)
    L1 = robot.L1;
    L2 = robot.L2;
    x = L1*cos(q1) + L2*cos(q1 + q2);
    y = L1*sin(q1) + L2*sin(q1 + q2);
end
 
function [q1, q2, reachable_geo] = IK(x, y, robot)
    L1 = robot.L1;
    L2 = robot.L2;
    c2 = (x^2 + y^2 - L1^2 - L2^2) / (2*L1*L2);
    if abs(c2) > 1
        % Point is outside the reachable annulus (|L1-L2| <= r <= L1+L2)
        reachable_geo = false;
        q1 = NaN;
        q2 = NaN;
        return;
    end
    reachable_geo = true;
    s2 = sqrt(1 - c2^2);       % elbow-up branch (s2 >= 0)
    q2 = atan2(s2, c2);
    q1 = atan2(y, x) - atan2(L2*s2, L1 + L2*c2);
end
 
function within = withinJointLimits(q1, q2, robot)
    within = (q1 >= robot.q1_limit(1)) && (q1 <= robot.q1_limit(2)) && ...
             (q2 >= robot.q2_limit(1)) && (q2 <= robot.q2_limit(2));
end
 
function [tau1_g_max, tau2_g_max] = gravityTorques(robot)
    g   = robot.g;
    m1  = robot.m1;  m2  = robot.m2;  mp = robot.mp;
    lc1 = robot.lc1; lc2 = robot.lc2;
    L1  = robot.L1;  L2  = robot.L2;
    tau2_g_max = g*(m2*lc2 + mp*L2);
    tau1_g_max = g*(m1*lc1 + m2*(L1 + lc2) + mp*(L1 + L2));
end
 
function [passFail, details] = driveFeasibility(tau_j_g_max, w_j_req, robot)
    ns    = robot.ns;
    G     = robot.G;
    eff_G = robot.eff_G;
    tau_j_req   = ns * tau_j_g_max;
    tau_m_req   = tau_j_req / (G * eff_G);
    w_m_req     = G * w_j_req;
    torque_ok = tau_m_req <= robot.t_m_rate;
    speed_ok  = w_m_req   <= robot.w_m_rate;
    passFail  = torque_ok && speed_ok;
    details.tau_j_req = tau_j_req;
    details.tau_m_req = tau_m_req;
    details.w_m_req   = w_m_req;
    details.torque_ok = torque_ok;
    details.speed_ok  = speed_ok;
end



