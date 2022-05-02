% AER1216 Fall 2021 
% Fixed Wing Project Code
%
% parameters.m
%
% Initialization file which generates and stores all required data into the 
% structure P, which is then stored in the workspace. Simulink model calls 
% on this function at the start of every simulation. Code structure adapted
% from Small Unmanned Aircraft: Theory and Practice by R.W. Beard and T. W. 
% McLain. 
% 
% Inputs: 
% N/A
%
% Outputs:
% P                 structure that contains all aerodynamic, geometric, and
%                   initial condition data for the aircraft and simulation.
%
% Last updated: Pravin Wedage 2021-11-09

%% TA NOTE
% An easy way to store parameters for use in simulink is through the use of
% a structure. For example, P.g = 9.81 stores the value of gravitational
% acceleration in the field g that is contained within the structure P.
% Anytime P is called anywhere in the simulation code, the value of P.g is
% accessible. 

%% Parameter Computation
% Initial Conditions
clear all
% compute trim conditions            
P.Va0 = 19;         % initial airspeed (also used as trim airspeed)
P.Va_trim = 15; 
P.Va = P.Va_trim;

%% CHANGES TO BE MADE

P.u_trim     = 18.236;
P.v_trim     = 0;
P.w_trim     = 5.01;
P.p_trim     = 0;
P.q_trim     = 0;
P.r_trim     = 0;
P.alpha_trim = 0.267; % ~15 degrees to radians
P.beta_trim  = 0;
P.dele_trim  = -0.25;
P.dela_trim  = 0;
P.delr_trim  = 0;
P.delt_trim  = 0.8311;
P.theta_trim = 0.267;
P.phi_trim   = 0;

P.Cprop      = 0.02;
P.kmotor     = 74; 

%%

P.gravity = 9.81;
P.g = 9.81; 
P.pho = 1.225;

% Aerosonde UAV Data
% physical parameters of airframe
P.m      = 13.5;    % kg
P.Ixx    = 0.8244;  % kg m^2
P.Iyy    = 1.135;   % kg m^2
P.Izz    = 1.759;   % kg m^2
P.Ixz    = 0.1204;  % kg m^2
P.S      = 0.55;    % m^2
P.b      = 2.8956;  % m
P.c      = 0.18944; % m  
P.Sprop  = 0.2027;  % m^2
P.e      = 0.9;     % oswald efficiency

% Longitudinal Aerodynamic Coefficients
P.CL_O     = 0.28;
P.CD_O     = 0.03;
P.Cm_O     = -0.02338;
P.CL_alpha = 3.45;
P.CD_alpha = 0.30;
P.Cm_alpha = -0.38;
P.CL_q     = 0;
P.CD_q     = 0;
P.Cm_q     = -3.6;
P.CL_dele  = -0.36;
P.CD_dele  = 0;
P.Cm_dele  = -0.5;
P.epsilon  = 0.1592;

% Lateral Aerodynamic Coefficients
P.CY_O     = 0;
P.Cl_O     = 0;
P.Cn_O     = 0;
P.CY_beta  = -0.98;
P.Cl_beta  = -0.12;
P.Cn_beta  = 0.25;
P.CY_p     = 0;
P.Cl_p     = -0.26;
P.Cn_p     = 0.022;
P.CY_r     = 0;
P.Cl_r     = 0.14;
P.Cn_r     = -0.35;
P.CY_dela  = 0;
P.Cl_dela  = 0.08;
P.Cn_dela  = 0.06;
P.CY_delr  = -0.17;
P.Cl_delr  = 0.105;
P.Cn_delr  = -0.032;

% Control Input limits 
P.delta_e_max = deg2rad(45); % assumed symmetric
P.delta_a_max = deg2rad(45); 
P.delta_r_max = deg2rad(25);

% Initial Conditions % connects with aircraft_dynamics.m, do not modify
% structure field names
P.pn0    = 0;  % initial North position
P.pe0    = 0;  % initial East position
P.pd0    = -2000;  % initial Down position (negative altitude)
P.u0     = P.Va0; % initial velocity along body x-axis
P.v0     = 0;  % initial velocity along body y-axis
P.w0     = 0;  % initial velocity along body z-axisu_trim
P.phi0   = 0;  % initial roll angle
P.theta0 = 0;  % initial pitch angle
P.psi0   = 0;  % initial yaw angle
P.p0     = 0;  % initial body frame roll rate
P.q0     = 0;  % initial body frame pitch rate
P.r0     = 0;  % initial body frame yaw rate
P.delta_e0 =0;
P.delta_a0 =0;
P.delta_r0 =0;
P.delta_t0 =0;
                        