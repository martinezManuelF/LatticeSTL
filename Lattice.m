%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lattice.m
%
% Manuel F. Martinez
%
% EE5322
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RESTORE STATE
close all;
clear all;
clc;

% UNITS
meters      = 1;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
degrees     = pi/180;

% OPEN FIGURE WINDOW
fig = figure('Color','w','DefaultFigureWindowStyle','docked');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LATTICE PARAMETERS
a = 1.0 * centimeters; % Lattice spacing
r = 0.3275 * centimeters; % Radius of spheres 

% GRID SIZE
Nx = 32;
Ny = Nx;
Nz = Nx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CREATE GRID, GENERATE LATTICE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GRID PARAMETERS
dx = a/Nx;
dy = a/Ny;
dz = a/Nz;

% GENERATE GRID AXES AND CENTER THEM
xa = [0 : Nx-1]*dx; xa = xa - mean(xa);
ya = [0 : Ny-1]*dy; ya = ya - mean(ya);
za = [0 : Nz-1]*dz; za = za - mean(za);

% GENERATE MESHGRID
[Y,X,Z] = meshgrid(ya,xa,za);

% GENERATE CORNER SPHERES
SP1 = [0.5 ; 0.5 ; 0.5] * a;
SP2 = [0.5 ; 0.5 ; -0.5] * a;
SP3 = [0.5 ; -0.5 ; 0.5] * a;
SP4 = [0.5 ; -0.5 ; -0.5] * a;
SP5 = [-0.5 ; 0.5 ; 0.5] * a;
SP6 = [-0.5 ; 0.5 ; -0.5] * a;
SP7 = [-0.5 ; -0.5 ; 0.5] * a;
SP8 = [-0.5 ; -0.5 ; -0.5] * a;
SP  = [SP1 SP2 SP3 SP4 SP5 SP6 SP7 SP8];

% GENERATE FACE SPHERES
SP1 = [0.5 ; 0 ; 0] * a;
SP2 = [-0.5 ; 0 ; 0] * a;
SP3 = [0 ; 0.5 ; 0] * a;
SP4 = [0 ; -0.5 ; 0] * a;
SP5 = [0 ; 0 ; 0.5] * a;
SP6 = [0 ; 0 ; -0.5] * a;
SP  = [SP SP1 SP2 SP3 SP4 SP5 SP6];

% GENERATE FCC SPHERES
SP1 = [-0.25 ; 0.25 ; 0.25] * a;
SP2 = [0.25 ; -0.25 ; 0.25] * a;
SP3 = [0.25 ; 0.25 ; -0.25] * a;
SP4 = [-0.25 ; -0.25 ; -0.25] * a;
SP  = [SP SP1 SP2 SP3 SP4];

% CONSTRUCT UNIT CELL
CELL = zeros(Nx,Ny,Nz); % Initialize unit cell
for sp = 1 : length(SP)
    offx    = SP(1,sp);     % Offset on X-Axis
    offy    = SP(2,sp);     % Offset on Y-Axis
    offz    = SP(3,sp);     % Offset on Z-Axis
    R       = (X + offx).^2 + (Y + offy).^2 + (Z + offz).^2;
    CELL    = CELL + (R < r^2);
end
CELL = ~CELL;

% GENERATE VERTICES AND FACES
[F,V]   = isosurface(xa,ya,za,CELL);
[F2,V2] = isocaps(xa,ya,za,CELL);
F = [F; F2+length(V(:,1))];
V = [V; V2];
patch('Faces',F,'Vertices',V,'LineStyle','none','FaceColor','b')
axis equal tight off;
camlight;
lighting phong;
view([170 15]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POST-PROCESS DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DISPLAY F AND V
display(['FACES     : ' num2str(length(F))]);
display(['VERTICES  : ' num2str(length(V))]);

% SAVE AS STL
stlwrite('DiamondLattice.stl',F,V,'mode','ascii');