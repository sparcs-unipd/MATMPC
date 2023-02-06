%------------------------------------------%

% high fidelity tilting quadrotor equations
% for MATMPC model generation

%------------------------------------------%


%% Dimensions

nx=21;  % No. of differential states (pos:3, quat:4, v:3, omega:3, alpha:4, w_bar2 (w_bar_square):4)
nu=8;   % No. of controls (dot_w_bar2:4, dot_alpha:4)
nz=0;   % No. of algebraic states
ny=28;  % No. of outputs (state + controls)
nyN=20; % No. of outputs at the terminal point (pos:3, quat_norm:1, quat_cos:1, v:3, omega:3, alpha:4, w_bar2:4)
np=4;   % No. of model parameters
nc=0;   % No. of general inequality constraints
ncN=0;  % No. of general inequality constraints
nbx = 8; % No. of bounds on states (alpha, w_bar2)
nbu = 8; % No. of bounds on controls (all)

% state and control bounds
nbx_idx = 14:21;  % indexs of states which are bounded (alpha)
nbu_idx = 1:8;  % indexs of controls which are bounded

%% create variables

import casadi.*

states   = SX.sym('states',nx,1);   % differential states
controls = SX.sym('controls',nu,1); % control input
alg      = SX.sym('alg',nz,1);      % algebraic states
params   = SX.sym('paras',np,1);    % parameters
refs     = SX.sym('refs',ny,1);     % references of the first N stages
refN     = SX.sym('refs',nyN,1);    % reference of the last stage
Q        = SX.sym('Q',ny,1);        % weighting matrix of the first N stages
QN       = SX.sym('QN',nyN,1);      % weighting matrix of the last stage
aux      = SX.sym('aux',ny,1);      % auxilary variable
auxN     = SX.sym('auxN',nyN,1);    % auxilary variable


%% Dynamics

g = 9.81; % gravity [m/s^2]

%%% Model Parameter loading

% look for existing parameter data
if ( exist('tiltq','var') ...
        && exist('body','var')  ...
        && exist('prop','var'))
    mR              = tiltq.m_b; % drone mass [Kg]
    JR              = body.I; % drone inertia [Kg*m^2]
    c_t_raw         = prop.k_t; % thrust coefficient [N/s^2]
    c_tau_raw       = prop.k_m; % drag coefficient [Nm/s^2]
    rotor_max_vel   = prop.w_max; % [rad/s]
    c_direction     = - prop.c_direction; % 1: CW , -1: CCW
    p_mot           = tiltq.p_mot;
else
    mR              = 2;
    JR              = 0.01*diag([1,1,2]);
    c_t_raw         = 1e-5;
    c_tau_raw       = 0.05*c_t_raw;
    rotor_max_vel   = 1000;
    c_direction     = [-1 -1 1 1];
    p_mot = [
         0.15 -0.15 0
        -0.15  0.15 0
         0.15  0.15 0
        -0.15 -0.15 0
        ];
end
theta           = atan2(p_mot(:,2),p_mot(:,1)); % arm direction angle
beta            = theta;
delta           = theta - beta;
l               = vecnorm(p_mot,2,2); % arm length

% normalized propeller spinning rate square
c_t = c_t_raw .* rotor_max_vel.^2;
c_tau = c_tau_raw .* rotor_max_vel.^2;

p           =states(1:3);
q           =states(4:7);
v           =states(8:10);
omega       =states(11:13);
alpha       =states(14:17);
w_bar2      =states(18:21); % propeller angular rate square
u           =controls(1:4); % propeller angular acceleration square
dot_alpha   =controls(5:8);
q_ref       =params(:);

% quaternion conversions
eta = q(1);
epsilon = q(2:4);
epsilon_x = [0, -epsilon(3), epsilon(2); epsilon(3), 0, -epsilon(1); -epsilon(2), epsilon(1), 0]; 
Rq = eye(3) + 2*eta*epsilon_x + 2*epsilon_x^2;  % rotation matrix, body to world
Mq = [eta, -epsilon'; epsilon, eta*eye(3)+epsilon_x]; % quaternion composition matrix


% thrust directions
F_temp = [...
    (sin(beta).*sin(alpha)).';...
    (-cos(beta).*sin(alpha)).';...
    (cos(alpha)).'];

% thrust matrix
F = 1/mR .* c_t .* F_temp;

% torque matrix
M_drag = c_tau .* F_temp * diag(c_direction);
M_thrust = c_t .* [...
    (sin(theta).*cos(alpha)).';...
    (-cos(theta).*cos(alpha)).';...
    (-cos(delta).*sin(alpha)).']...
    *diag(l);
M = M_drag + M_thrust;

% explicit ODE RHS
x_dot = [
    v
    1/2*Mq*[0;omega]
    Rq*F*w_bar2-[0;0;g]
    JR \ M*w_bar2
    dot_alpha
    u
    ];

% algebraic function
z_fun = [];

% implicit ODE: impl_f = 0
xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;
     
%% Objectives and constraints

q_err = Mq*[q_ref(1); -q_ref(2:4)];


% inner objectives
h = [
    p
    q_err(2:end)
    v
    omega
    alpha
    w_bar2
    u
    dot_alpha
    ];
hN = h(1:nyN);

% outer objectives
obji = 0.5*(h-refs)'*diag(Q)*(h-refs);
objN = 0.5*(hN-refN)'*diag(QN)*(hN-refN);

obji_GGN = 0.5*(aux-refs)'*(aux-refs);
objN_GGN = 0.5*(auxN-refN)'*(auxN-refN);

% general inequality path constraints
general_con = [];
general_con_N = [];

%% NMPC discretizing time length [s]

Ts_st = 0.01; % shooting interval time

%% Model Parameters exportation
modelParameters = struct(...
    'name', 'highFidelityTiltingQuadrotor', ...
    'g', g, ...
    'mR', mR, ...
    'JR', JR, ...
    'c_t_raw', c_t_raw, ...
    'c_tau_raw', c_tau_raw, ...
    'rotor_max_vel', rotor_max_vel, ...
    'c_direction', c_direction, ...
    'p_mot', p_mot, ...
    'theta', theta,...
    'beta', beta, ...
    'delta', delta, ...
    'l', l);
save('modelParameters.mat', 'modelParameters');
