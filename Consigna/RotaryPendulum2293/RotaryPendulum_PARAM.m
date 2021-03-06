%RotaryPendulum_PARAM    Parameter file to initialize the rotary pendulum
% -------------------------------------------------------------------------
% inputs    : -
% outputs   : -
% -------------------------------------------------------------------------
% Copyright 2015 The MathWorks, Inc.
% -------------------------------------------------------------------------

%% ------------------------------------------------------------------------
%  model environment properties
%  ------------------------------------------------------------------------
g = 9.81;

%% ------------------------------------------------------------------------
%  Base (just for PhysMod geometry)
%  ------------------------------------------------------------------------
radius_base= 40e-3;
length_base= 10e-3;

%% ------------------------------------------------------------------------
%  Pin (just for PhysMod geometry)
%  ------------------------------------------------------------------------
radius_vert_cylinder = 10e-3;
length_vert_cylinder = 200e-3;

%% ------------------------------------------------------------------------
%  RotaryArm
%  ------------------------------------------------------------------------
M = 20.3e-3;
R = 109e-3;

me    = 0.012;              % mass encoder
ma    = M - me;

% PhysMod geometry
horizontal_arm_geo = [109e-3 radius_vert_cylinder*2 radius_vert_cylinder*2];

% Tradiotioanl approach parameters
re    = 0.02;               % encoder offset
Ino_e = ma*R^2/3;           % corrected I without encoder
Itot  = Ino_e + me*re^2;    % corrected I with encoder
I     = Itot;

%% ------------------------------------------------------------------------
%  Pendulum
%  ------------------------------------------------------------------------
m = 3.3e-3;
l = 183.2e-3;

% PhysMod geometry
pendulum_geo = [5e-3 20e-3 l];

%% ------------------------------------------------------------------------
%  Setup Pendulum_CAD_Control initial values
%  ------------------------------------------------------------------------
friction_on_off = 1;
controller_on_off = 0;

gain_motor = 5.67;
