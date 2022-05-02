function aircraft_dynamics(block)
%MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the 
%   name of your S-function.
%
%   It should be noted that the MATLAB S-function is very similar
%   to Level-2 C-Mex S-functions. You should be able to get more
%   information for each of the block methods by referring to the
%   documentation for C-Mex S-functions.
%
%   Copyright 2003-2010 The MathWorks, Inc.

% AER1216 Fall 2021 
% Fixed Wing Project Code
%
% aircraft_dynamics.m
%
% Fixed wing simulation model file, based on the Aerosonde UAV, with code
% structure adapted from Small Unmanned Aircraft: Theory and Practice by 
% R.W. Beard and T. W. McLain. 
% 
% Inputs: 
% delta_e           elevator deflection [deg]
% delta_a           aileron deflection [deg]
% delta_r           rudder deflection [deg]
% delta_t           normalized thrust []
%
% Outputs:
% pn                inertial frame x (north) position [m]
% pe                inertial frame y (east) position [m]
% pd                inertial frame z (down) position [m]
% u                 body frame x velocity [m/s]
% v                 body frame y velocity [m/s]
% w                 body frame z velocity [m/s]
% phi               roll angle [rad]
% theta             pitch angle [rad]
% psi               yaw angle [rad]
% p                 roll rate [rad/s]
% q                 pitch rate [rad/s]
% r                 yaw rate [rad/s]
%
% Last updated: Pravin Wedage 2021-11-09

%% TA NOTE
% The code segements you must modify are located in the derivatives
% function in this .m file. Modify other sections at your own risk. 


%
% The setup method is used to set up the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.
%
setup(block);

end 


%% Function: setup ===================================================
% Abstract:
%   Set up the basic characteristics of the S-function block such as:
%   - Input ports
%   - Output ports
%   - Dialog parameters
%   - Options
%
%   Required         : Yes
%   C-Mex counterpart: mdlInitializeSizes
%
function setup(block)

% Register number of ports
block.NumInputPorts  = 1;
block.NumOutputPorts = 1;

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
for i = 1:block.NumInputPorts
    block.InputPort(i).Dimensions        = 4;
    block.InputPort(i).DatatypeID  = 0;  % double
    block.InputPort(i).Complexity  = 'Real';
    block.InputPort(i).DirectFeedthrough = false; % important to be false 
end

% Override output port properties
for i = 1:block.NumOutputPorts
    block.OutputPort(i).Dimensions       = 12;
    block.OutputPort(i).DatatypeID  = 0; % double
    block.OutputPort(i).Complexity  = 'Real';
%     block.OutputPort(i).SamplingMode = 'Sample';
end

% Register parameters
block.NumDialogPrms     = 1;
P = block.DialogPrm(1).Data; % must duplicate this line in each function

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [0 0];

% Register multiple instances allowable
% block.SupportMultipleExecInstances = true;

% Register number of continuous states
block.NumContStates = 12;

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

% -----------------------------------------------------------------
% The MATLAB S-function uses an internal registry for all
% block methods. You should register all relevant methods
% (optional and required) as illustrated below. You may choose
% any suitable name for the methods and implement these methods
% as local functions within the same file. See comments
% provided for each function for more information.
% -----------------------------------------------------------------

% block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup); % discrete states only
block.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);
block.RegBlockMethod('InitializeConditions',    @InitializeConditions);
% block.RegBlockMethod('Start',                   @Start); % Initialize Conditions is used
block.RegBlockMethod('Outputs',                 @Outputs); % Required
% block.RegBlockMethod('Update',                  @Update); % only required for discrete states
block.RegBlockMethod('Derivatives',             @Derivatives); % Required for continuous states
block.RegBlockMethod('Terminate',               @Terminate); % Required

end 


%% PostPropagationSetup:
%   Functionality    : Setup work areas and state variables. Can
%                      also register run-time methods here
%   Required         : No
%   C-Mex counterpart: mdlSetWorkWidths
%
function DoPostPropSetup(block)
block.NumDworks = 1;
  
  block.Dwork(1).Name            = 'x1';
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 0;      % double
  block.Dwork(1).Complexity      = 'Real'; % real
  block.Dwork(1).UsedAsDiscState = true;

end


%% InitializeConditions:
%   Functionality    : Called at the start of simulation and if it is 
%                      present in an enabled subsystem configured to reset 
%                      states, it will be called when the enabled subsystem
%                      restarts execution to reset the states.
%   Required         : No
%   C-MEX counterpart: mdlInitializeConditions
%
function InitializeConditions(block)

% Rename parameters
P = block.DialogPrm(1).Data; % must duplicate this line in each function

% Initialize continuous states
block.ContStates.Data(1) = P.pn0; 
block.ContStates.Data(2) = P.pe0;
block.ContStates.Data(3) = P.pd0;
block.ContStates.Data(4) = P.u0;
block.ContStates.Data(5) = P.v0;
block.ContStates.Data(6) = P.w0;
block.ContStates.Data(7) = P.phi0;
block.ContStates.Data(8) = P.theta0;
block.ContStates.Data(9) = P.psi0;
block.ContStates.Data(10) = P.p0;
block.ContStates.Data(11) = P.q0;
block.ContStates.Data(12) = P.r0;

end 

%% Start:
%   Functionality    : Called once at start of model execution. If you
%                      have states that should be initialized once, this 
%                      is the place to do it.
%   Required         : No
%   C-MEX counterpart: mdlStart
%
function Start(block)

block.Dwork(1).Data = 0;

end 

%% Input Port Sampling Method:
function SetInpPortFrameData(block, idx, fd)
  
  block.InputPort(idx).SamplingMode = 'Sample';
  for i = 1:block.NumOutputPorts
    block.OutputPort(i).SamplingMode  = 'Sample';   
  end
end

%% Outputs:
%   Functionality    : Called to generate block outputs in
%                      simulation step
%   Required         : Yes
%   C-MEX counterpart: mdlOutputs
%
function Outputs(block)

temp_mat = zeros(block.NumContStates,1); % thirteen states
for i = 1:block.NumContStates
     temp_mat(i) = block.ContStates.Data(i);
end

block.OutputPort(1).Data = temp_mat; % states

% for i = 1:block.NumOutputPorts
%     block.OutputPort(1).Data(i) = block.ContStates.Data(i);
% end

end 


%% Update:
%   Functionality    : Called to update discrete states
%                      during simulation step
%   Required         : No
%   C-MEX counterpart: mdlUpdate
%
function Update(block)

block.Dwork(1).Data = block.InputPort(1).Data;

end 


%% Derivatives:
%   Functionality    : Called to update derivatives of
%                      continuous states during simulation step
%   Required         : No
%   C-MEX counterpart: mdlDerivatives
%
function Derivatives(block)

% Rename parameters
P = block.DialogPrm(1).Data; % must duplicate this line in each function

% compute inertial constants
T     = P.Ixx*P.Izz-P.Ixz^2;
T1    = (P.Ixz*(P.Ixx-P.Iyy+P.Izz))/T;
T2    = (P.Izz*(P.Izz-P.Iyy)+P.Ixz^2)/T;
T3    = P.Izz/T;
T4    = P.Ixz/T;
T5    = (P.Izz-P.Ixx)/P.Iyy;
T6    = P.Ixz/P.Iyy;
T7    = (P.Ixx*(P.Ixx-P.Iyy)+P.Ixz^2)/T;
T8    = P.Ixx/T;

% map states and inputs
pn    = block.ContStates.Data(1);
pe    = block.ContStates.Data(2);
pd    = block.ContStates.Data(3);
u     = block.ContStates.Data(4);
v     = block.ContStates.Data(5);
w     = block.ContStates.Data(6);
phi   = block.ContStates.Data(7);
theta = block.ContStates.Data(8);
psi   = block.ContStates.Data(9);
p     = block.ContStates.Data(10);
q     = block.ContStates.Data(11);
r     = block.ContStates.Data(12);
delta_e = block.InputPort(1).Data(1)*pi/180 ; % converted inputs to radians
delta_a = block.InputPort(1).Data(2)*pi/180 ; % converted inputs to radians
delta_r = block.InputPort(1).Data(3)*pi/180 ; % converted inputs to radians
delta_t = block.InputPort(1).Data(4);


% Air Data 
Va    = sqrt(u^2 + v^2 + w^2);
alpha = atan(w/u);
beta  = atan(v/sqrt(u^2 + w^2));

% rotation matrix
rot_xyz_dot    = [cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);
                  cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);
                 -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)];


% IMP: pndot(xdot) and pedot(ydot) are calculated using non-linear dynamics             
temp  = rot_xyz_dot*[u; v; w];
pndot = temp(1);
pedot = temp(2);


%% Longitudinal System

% State-Space Model Coefficients Table 5.2 (Page 86, Chapter 5)
CX_O     = -cos(alpha)*P.CD_O + sin(alpha)*P.CL_O;
CX_alpha = -cos(alpha)*P.CD_alpha + sin(alpha)*P.CL_alpha;
CX_dele  = -cos(alpha)*P.CD_dele + sin(alpha)*P.CL_dele;
CX_q     = -cos(alpha)*P.CD_q + sin(alpha)*P.CL_q;
CZ_O     = -cos(alpha)*P.CL_O - sin(alpha)*P.CD_O;
CZ_alpha = -cos(alpha)*P.CL_alpha - sin(alpha)*P.CD_alpha;
CZ_dele  = -cos(alpha)*P.CL_dele - sin(alpha)*P.CD_dele;
CZ_q     = -cos(alpha)*P.CL_q - sin(alpha)*P.CD_q;

Xu       = (P.u_trim*P.pho*P.S/P.m)*(CX_O + CX_alpha*P.alpha_trim + CX_dele*P.dele_trim) ...
            -(P.pho*P.S*P.w_trim*CX_alpha/2*P.m) + ...
             (P.pho*P.S*P.c*CX_q*P.u_trim*P.q_trim/4*P.m*P.Va_trim) - ...
             (P.pho*P.Sprop*P.Cprop*P.u_trim/P.m);

Xw       = -P.q_trim + ...
           (P.w_trim*P.pho*P.S/P.m)*(CX_O + CX_alpha*P.alpha_trim + CX_dele*P.dele_trim) + ...
           (P.pho*P.S*P.u_trim*CX_alpha/2*P.m) + ...
           (P.pho*P.S*P.c*CX_q*P.w_trim*P.q_trim)/(4*P.m*P.Va_trim) - ...
           (P.pho*P.Sprop*P.Cprop*P.w_trim/P.m);
           
Xq       = -P.w_trim + (P.pho*P.Va_trim*P.S*P.c*CX_q/4*P.m);

X_dele   = P.pho*(P.Va_trim^2)*P.S*CX_dele/2*P.m;

X_delt   = P.pho*P.Sprop*P.Cprop*(P.kmotor^2)*P.delt_trim/P.m;

Zu       = P.q_trim + ...
          (P.u_trim*P.pho*P.S/P.m)*(CZ_O + CZ_alpha*P.alpha_trim + CZ_dele*P.dele_trim) - ...
          (P.pho*P.S*P.w_trim*CZ_alpha/2*P.m) + ...
          (P.pho*P.S*P.c*CZ_q*P.u_trim*P.q_trim)/(4*P.m*P.Va_trim);
          
Zw       = (P.w_trim*P.pho*P.S/P.m)*(CZ_O + CZ_alpha*P.alpha_trim + CZ_dele*P.dele_trim) + ...
           (P.pho*P.S*P.u_trim*CZ_alpha/2*P.m) + ...
           (P.pho*P.S*P.c*CZ_q*P.w_trim*P.q_trim)/(4*P.m*P.Va_trim);
       
Zq       = P.u_trim + ...
          (P.pho*P.S*P.c*CZ_q*P.Va_trim)/(4*P.m);
      
Z_dele   = (P.pho*P.S*CZ_dele*(P.Va_trim^2))/(2*P.m);

Mu       = (P.u_trim*P.pho*P.S*P.c/P.Iyy)*(P.Cm_O + P.Cm_alpha*P.alpha_trim + P.Cm_dele*P.dele_trim) - ...
           (P.pho*P.S*P.c*P.w_trim*P.Cm_alpha/2*P.Iyy) + ...
           (P.pho*P.S*(P.c^2)*P.Cm_q*P.u_trim*P.q_trim)/(4*P.Iyy*P.Va_trim);
       
Mw       = (P.w_trim*P.pho*P.S*P.c/P.Iyy)*(P.Cm_O + P.Cm_alpha*P.alpha_trim + P.Cm_dele*P.dele_trim) + ...
           (P.pho*P.S*P.c*P.u_trim*P.Cm_alpha/2*P.Iyy) + ...
           (P.pho*P.S*(P.c^2)*P.Cm_q*P.w_trim*P.q_trim)/(4*P.Iyy*P.Va_trim);
       
Mq       = (P.pho*P.S*(P.c^2)*P.Cm_q*P.Va_trim)/(4*P.Iyy);

M_dele   = (P.pho*P.S*P.c*P.Cm_dele*(P.Va_trim^2))/(2*P.Iyy);

% Longitudinal Linear Equations
A        = [Xu Xw Xq -P.g*cos(P.theta_trim) 0;
            Zu Zw Zq -P.g*sin(P.theta_trim) 0;
            Mu Mw Mq 0 0;
            0 0 1 0 0;
            sin(P.theta_trim) -cos(P.theta_trim) 0 P.u_trim*cos(P.theta_trim)+P.w_trim*sin(P.theta_trim) 0];

h        = -pd;
X        = [u; w; q; theta; h];

B        = [X_dele X_delt;
            Z_dele 0;
            M_dele 0;
            0 0;
            0 0];
        
U        = [delta_e; delta_t];

Xdot     = A*X + B*U;
udot     = Xdot(1);
wdot     = Xdot(2);
qdot     = Xdot(3);
thetadot = Xdot(4);
pddot    = -Xdot(5);   % because pd = -h

%% Lateral System

% Coefficient Definitions
Cp_O     = T3*P.Cl_O + T4*P.Cn_O;
Cp_beta  = T3*P.Cl_beta + T4*P.Cn_beta;
Cp_p     = T3*P.Cl_p + T4*P.Cn_p;
Cp_r     = T3*P.Cl_r + T4*P.Cn_r;
Cp_dela  = T3*P.Cl_dela + T4*P.Cn_dela;
Cp_delr  = T3*P.Cl_delr + T4*P.Cn_delr;
Cr_O     = T4*P.Cl_O + T8*P.Cn_O;
Cr_beta  = T4*P.Cl_beta + T8*P.Cn_beta;
Cr_p     = T4*P.Cl_p + T8*P.Cn_p;
Cr_r     = T4*P.Cl_r + T8*P.Cn_r;
Cr_dela  = T4*P.Cl_dela + T8*P.Cn_dela;
Cr_delr  = T4*P.Cl_delr + T8*P.Cn_delr;

% Lateral State-Space Model Coefficients Table 5.1 (Page 81, Chapter 5)
Yv       = ((P.pho*P.S*P.b*P.v_trim)/(4*P.m*P.Va_trim))*(P.CY_p*P.p_trim + P.CY_r*P.r_trim) + ...
           (P.pho*P.S*P.v_trim/P.m)*(P.CY_O + P.CY_beta*P.beta_trim + P.CY_dela*P.dela_trim + P.CY_delr*P.delr_trim) + ...
           (P.pho*P.S*P.CY_beta/2*P.m)*sqrt(P.u_trim^2 + P.w_trim^2);


Yp       = P.w_trim + (P.pho*P.Va_trim*P.S*P.b/(4*P.m))*P.CY_p;

Yr       = -P.u_trim + (P.pho*P.Va_trim*P.S*P.b/(4*P.m))*P.CY_r;

Y_dela   = (P.pho*(P.Va_trim^2)*P.S/(2*P.m))*P.CY_dela;

Y_delr   = (P.pho*(P.Va_trim^2)*P.S/(2*P.m))*P.CY_delr;

Lv       = ((P.pho*P.S*(P.b^2)*P.v_trim)/(4*P.Va_trim))*(Cp_p*P.p_trim + Cp_r*P.r_trim) + ...
           (P.pho*P.S*P.b*P.v_trim)*(Cp_O + Cp_beta*P.beta_trim + Cp_dela*P.dela_trim + Cp_delr*P.delr_trim) + ...
           (P.pho*P.S*P.b*Cp_beta/2)*sqrt(P.u_trim^2 + P.w_trim^2);

Lp       = T1*P.q_trim + (P.pho*P.Va_trim*P.S*(P.b^2)/4)*Cp_p;

Lr       = -T2*P.q_trim + (P.pho*P.Va_trim*P.S*(P.b^2)/4)*Cp_r;

L_dela   = (P.pho*(P.Va_trim^2)*P.S*P.b/2)*Cp_dela;

L_delr   = (P.pho*(P.Va_trim^2)*P.S*P.b/2)*Cp_delr;

Nv       = ((P.pho*P.S*(P.b^2)*P.v_trim)/(4*P.Va_trim))*(Cr_p*P.p_trim + Cr_r*P.r_trim) + ...
           (P.pho*P.S*P.b*P.v_trim)*(Cr_O + Cr_beta*P.beta_trim + Cr_dela*P.dela_trim + Cr_delr*P.delr_trim) + ...
           (P.pho*P.S*P.b*Cr_beta/2)*sqrt(P.u_trim^2 + P.w_trim^2);
       
Np       = T7*P.q_trim + (P.pho*P.Va_trim*P.S*(P.b^2)/4)*Cr_p;

Nr       = -T1*P.q_trim + (P.pho*P.Va_trim*P.S*(P.b^2)/4)*Cr_r;

N_dela   = (P.pho*(P.Va_trim^2)*P.S*P.b/2)*Cr_dela;

N_delr   = (P.pho*(P.Va_trim^2)*P.S*P.b/2)*Cr_delr;
       
% Lateral Linear Equations
A        = [Yv Yp Yr P.g*cos(P.theta_trim)*cos(P.phi_trim) 0;
            Lv Lp Lr 0 0;
            Nv Np Nr 0 0;
            0 1 cos(P.phi_trim)*tan(P.theta_trim) P.q_trim*cos(P.phi_trim)*tan(P.theta_trim)-P.r_trim*sin(P.phi_trim)*tan(P.theta_trim) 0;
            0 0 cos(P.phi_trim)*sec(P.theta_trim) P.p_trim*cos(P.phi_trim)*sec(P.theta_trim)-P.r_trim*sin(P.phi_trim)*sec(P.theta_trim) 0];

X        = [v; p; r; phi; psi];

B        = [Y_dela Y_delr;
            L_dela L_delr;
            N_dela N_delr;
            0 0;
            0 0];
        
U        = [delta_a; delta_r];

Xdot     = A*X + B*U;
vdot     = Xdot(1);
pdot     = Xdot(2);
rdot     = Xdot(3);
phidot   = Xdot(4);
psidot   = Xdot(5);   

%%

% map derivatives
block.Derivatives.Data(1) = pndot;
block.Derivatives.Data(2) = pedot;
block.Derivatives.Data(3) = pddot;
block.Derivatives.Data(4) = udot;
block.Derivatives.Data(5) = vdot;
block.Derivatives.Data(6) = wdot;
block.Derivatives.Data(7) = phidot;
block.Derivatives.Data(8) = thetadot;
block.Derivatives.Data(9) = psidot;
block.Derivatives.Data(10)= pdot;
block.Derivatives.Data(11)= qdot;
block.Derivatives.Data(12)= rdot;

end 


%% Terminate:
%   Functionality    : Called at the end of simulation for cleanup
%   Required         : Yes
%   C-MEX counterpart: mdlTerminate
%
function Terminate(block)

end 


