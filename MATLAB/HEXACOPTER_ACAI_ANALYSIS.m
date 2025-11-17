%% HEXACOPTER ACAI ANALYSIS
% This script computes the Available Control Authority Index (ACAI)
% for a simple hexacopter (PNPNPN) and then sweeps over degraded thrust
% levels for each rotor.
%
% Author: Zahra Motahari Rad
% Description:
%   1) Define hexacopter parameters and control effectiveness matrix.
%   2) Compute baseline controllability margin rho_max.
%   3) Sweep over scaled thrust levels (0.0, 0.5, 0.8, 0.9, 1.0) for each rotor.
%   4) Compute the normalized margin rhobar for each thrust combination.
%
% Notation:
%   T  : total thrust
%   L  : roll moment
%   M  : pitch moment
%   N  : yaw moment
%
% This script is intended as a simple reproducible example to accompany
% controllability analysis for the PNPNPN hexacopter layout.

clc;
clear;

%% Design parameters
param.m_a   = 3.87;     % total mass [kg]
param.d     = 0.40;     % rotor arrangement radius [m]
param.f_max = 14.5;     % maximum rotor thrust [N]
param.k     = 0.03;     % reactive torque to thrust ratio [m]

%% Disturbance vector (weight compensation)
% Disturbance in [T L M N]' that must be balanced by the rotors
G  = -[-param.m_a * 9.79, 0, 0, 0]';  % [T L M N]'
u0 = G;                                % Disturbance -> compensation required

%% Control effectiveness matrix for PNPNPN hexacopter
% Rotor sequence PNPNPN for yaw direction
th =           [1,  1,  1,  1,  1,  1];                  % thrust
b1 = sqrt(3)/2 * [0, -1, -1,  0,  1,  1];                % roll
b2 =           [1,  0.5, -0.5, -1, -0.5, 0.5];           % pitch
b3 =           [1, -1,  1, -1,  1, -1];                  % yaw

disp('Simple hexacopter (PNPNPN)');

% Control effectiveness matrix B_f (mapping rotor forces to [T L M N]')
param.B_f = diag([1, param.d, param.d, param.k]) * [th; b1; b2; b3];

% Effector limits (all rotors healthy, full authority)
u_min = zeros(size(param.B_f, 2), 1);              % lower thrust bound [0 ... 0]'
u_max = param.f_max * ones(size(param.B_f, 2), 1); % upper thrust bound [f_max ... f_max]'

%% Helper: compute ACAI margin for given bounds
compute_margin = @(u_min_local, u_max_local) ...
    local_acai_margin(param.B_f, u_min_local, u_max_local, G);

%% 1) Baseline controllability margin (all rotors healthy)
rho_max = compute_margin(u_min, u_max);

%% 2) Thrust degradation sweep
% Thrust scaling levels for each rotor
th_levels = [0.0, 0.5, 0.8, 0.9, 1.0];

% Number of combinations (5^6 = 15625)
n_levels  = numel(th_levels);
n_rotors  = 6;
n_cases   = n_levels^n_rotors;

% Preallocate result arrays
rhobar   = zeros(n_cases, 1);
Thrust1  = zeros(n_cases, 1);
Thrust2  = zeros(n_cases, 1);
Thrust3  = zeros(n_cases, 1);
Thrust4  = zeros(n_cases, 1);
Thrust5  = zeros(n_cases, 1);
Thrust6  = zeros(n_cases, 1);

counter = 1;

for T1 = th_levels
    for T2 = th_levels
        for T3 = th_levels
            for T4 = th_levels
                for T5 = th_levels
                    for T6 = th_levels

                        % Current thrust scaling vector for each rotor
                        th_scale = [T1, T2, T3, T4, T5, T6]';

                        % Update effector limits for this case
                        u_min_case = zeros(size(param.B_f, 2), 1);
                        u_max_case = param.f_max * th_scale;

                        % Compute margin for this degraded case
                        rho = compute_margin(u_min_case, u_max_case);

                        % Normalized margin with respect to healthy case
                        rhobar(counter, 1) = rho / rho_max;

                        % Save thrust scaling levels
                        Thrust1(counter, 1) = T1;
                        Thrust2(counter, 1) = T2;
                        Thrust3(counter, 1) = T3;
                        Thrust4(counter, 1) = T4;
                        Thrust5(counter, 1) = T5;
                        Thrust6(counter, 1) = T6;

                        counter = counter + 1;

                    end
                end
            end
        end
    end
end

%% Optional: pack results in a struct for saving
results.rhobar  = rhobar;
results.T1      = Thrust1;
results.T2      = Thrust2;
results.T3      = Thrust3;
results.T4      = Thrust4;
results.T5      = Thrust5;
results.T6      = Thrust6;
results.rho_max = rho_max;

% Uncomment to save:
% save('hex_paca_results.mat', 'results');

disp('Sweep finished.');

%% ------------------------------------------------------------------------
%  Local function: ACAI margin computation
% -------------------------------------------------------------------------
function rho = local_acai_margin(B_f, u_min, u_max, G)
    % Compute ACAI  margin for a given control effectiveness matrix B_f,
    % effector bounds u_min/u_max, and disturbance G.

    delta = u_max - u_min;       % linear effector contribution magnitude
    fc    = (u_max + u_min) / 2; % center of envelope in effector space
    Fc    = B_f * fc;            % center of envelope in force/moment space

    [n_Bf, m_Bf] = size(B_f);
    M  = 1:m_Bf;
    S1 = nchoosek(M, n_Bf - 1);  % combinations of (n-1) effectors

    n_segments = size(S1, 1);
    d          = zeros(n_segments, 1);

    for i = 1:n_segments
        choose = S1(i, :);       % combination of (n-1) effectors

        % B1: hypersegment row-space
        B1 = B_f(:, choose);

        % xi: unit vector orthogonal to hypersegment column-space
        xi = null(B1');
        xi = xi(:, 1) / norm(xi(:, 1));

        % Remaining effectors and their ranges
        B2f    = B_f;
        B2f(:, choose) = [];
        delta2        = delta;
        delta2(choose) = [];

        % Sum of projected absolute contribution
        l = abs(xi' * B2f) * abs(delta2);

        % Absolute projected distance from G to Fc
        g = abs(xi' * (Fc - G));

        % Difference between G-to-Fc and hyperplane segment-to-Fc
        d(i) = l / 2 - g;
    end

    % Margin: minimum distance over all hyperplane segments
    rho = min(d);
end

