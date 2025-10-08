%% 3D Conservative Stokes / Marker in cell / plastic iterations
% VisoElastoPlasic approach
% by Jonas, March/April 2012
clear all; 

addpath('SOURCE_FILES')                        % FUNCTIONS
mkdir OUTPUT

create_output       = 100;  % producing output files every x timestep
create_breakpoint   = 200;  % producing breakpoint files every x timestep

brutus              = 1;    % 1 if you run it on cluster (without plotting every timestep)
Temperature         = 0;
Powerlaw            = 0;
Diffusion_creep     = 0;
Surface_diffusion   = 1;
Beam_function       = 0;    % 1: beam function of elastic bending at the bottom. Requires material below the bottom and no boundary velocities
delta_prograd       = 0;    % 1: delta progradation on the left side of the model
Variable_grid       = 0;
variable_lambda     = 0;

if brutus == 0
    clf;
    load Colormaps;
end

% Length of model
Lx   = 100e3;
Ly   =  15e3;

% number of nodes
if Variable_grid == 1
    % dx = [2000*ones(1,140) 2000:-100:500 500*ones(1,800) 500:100:2000 2000*ones(1,240)];
    % dy = [500*ones(1,140) 500:100:2000 2000*ones(1,60)];
    % for i = 1:length(dx)
    %     x(i+1)  =  sum(dx(1:i));
    % end
    % for i = 1:length(dy)
    %     y(i+1)  =  sum(dy(1:i));
    % end
    k=5000;
    for i = 2:22
        k(i) = k(i-1).*0.9;
    end
    f=sum(k(2:end));
    f=f/40e3;
    xxx = k(2:end)/f;
    xxxx = flip(xxx);
    for i = 1:length(xxx)
        xxx1(i) = sum(xxx(1:i));
        xxx2(i) = sum(xxxx(1:i));
    end
    x = [0:5e3:500e3 500e3+xxx1 540.5e3:500:960e3 960e3+xxx2 1005e3:5e3:1500e3];
    y = [0:500:70e3 70e3+xxx2 115e3:5e3:250e3];

    dx = x(2:end)-x(1:end-1);
    dy = y(2:end)-y(1:end-1);

    nx      =   length(x);
    ny      =   length(y);
    p_node  =   ceil(nx/2);
else
    nx      =   1001;
    ny      =    151;
    p_node  =    ceil(nx/2);   % between 3 and nx-2
end
clear xxx xxx1 xxx2 xxxx i k f;

% markers per node
mx  =   5;
my  =   5;

        % Phase    Visc       EModul    Dens     Phi     PhiW      Coh         CohW      lamb      Hr       Ta      Tb       cp      Temp       kk       n       Q      DifCreep    grains      m      Qd
ROCKS= [    1      1e18        1e11        1       0       1       1e20        1e20        0        0        0       0      3e6       273      200       1       0             0        0       0       0   ;   % AIR
            2      1e18        1e11      1e3       0       1       1e20        1e20        0        0        0       0      3e3       273      200       1       0             0        0       0       0   ;   % WATER
%             3   1.9953e21      1e12     2500      30      30       10e6         1e6      0.90    2e-6    2e-5   45e-13      1e3       273      2.5      3.0  220e3      1.585e18    50e-6      -2   220e3   ;   % Sediment 1; Quartz (Brodie and Rutter, 2000)
            3   1.5849e35      1e11     2700      30      25        5e6         1e6      0.9    2e-6    2.5e-5   45e-13      1e3       273      2.5      4.0  135e3      1.585e18     1e-3      -2   220e3   ;   % Sediment 1; Quartz (Hirth et al. 2001; Brodie and Rutter, 2000)
            4   1.5849e35      1e11     2700      30      25        5e6         1e6      0.4     2e-6    2e-5   45e-13      1e3       273      2.5      4.0  135e3      1.585e18     1e-3      -2   220e3   ;   % Sediment 1; Quartz (Hirth et al. 2001; Brodie and Rutter, 2000)
            5   1.5849e35      1e11     2700      30      25        5e6         1e6      0.4     2e-6    2e-5   45e-13      1e3       273      2.5      4.0  135e3      1.585e18     1e-3      -2   220e3   ;   % Sediment 1; Quartz (Hirth et al. 2001; Brodie and Rutter, 2000)
            6   1.5849e35      1e11     2700      30      25        5e6         1e6      0.4     2e-6    2e-5   45e-13      1e3       273      2.5      4.0  135e3      1.585e18     1e-3      -2   220e3   ;   % Sediment 1; Quartz (Hirth et al. 2001; Brodie and Rutter, 2000)
            7   1.5849e35      1e11     2700      30      25        5e6         1e6      0.4     2e-6    2e-5   45e-13      1e3       273      2.5      4.0  135e3      1.585e18     1e-3      -2   220e3   ;   % Sediment 1; Quartz (Hirth et al. 2001; Brodie and Rutter, 2000)
            8   1.5849e35      1e11     2700      30      25        5e6         1e6      0.4     2e-6    2e-5   45e-13      1e3       273      2.5      4.0  135e3      1.585e18     1e-3      -2   220e3   ;   % Sediment 1; Quartz (Hirth et al. 2001; Brodie and Rutter, 2000)
            9   1.5849e35      1e11     2800      30      20       10e6         1e6      0.4     2e-6    2e-5   45e-13      1e3       273      2.0      4.0  135e3      1.585e18     1e-2      -2   220e3   ;   % Upper Crust; Quartz (Hirth et al. 2001; Brodie and Rutter, 2000)
           10   1.5849e35      1e11     2800      30      20       10e6         1e6      0.4     2e-6    2e-5   45e-13      1e3       273      2.0      4.0  135e3      1.585e18     1e-2      -2   220e3   ;   % Upper Crust; Quartz (Hirth et al. 2001; Brodie and Rutter, 2000)
%            11   1.9952e11      1e11     2800      30      15       25e6         1e6      0.0     2e-7    2e-5   45e-13      1e3       273      1.6      3.0  648e3     7.9428e11     1e-2      -3   467e3   ;   % Lower Crust; Anorthite (Rybacki and Dresen, 2000)
%            12   1.9952e11      1e11     2800      30      15       25e6         1e6      0.0     2e-7    2e-5   45e-13      1e3       273      1.6      3.0  648e3     7.9428e11     1e-2      -3   467e3   ;   % Lower Crust; Anorthite (Rybacki and Dresen, 2000)
%            11   3.1623e19      1e11     2700      35      35       25e6        25e6      0.0     2e-7    2e-5   45e-13      1e3       273      1.6      3.0  235e3     7.9433e22     2e-2      -3   153e3   ;   % Wet Anorthite60 (Rybacki and Dresen, 2004)
%            12   3.1623e19      1e11     2700      35      35       25e6        25e6      0.0     2e-7    2e-5   45e-13      1e3       273      1.6      3.0  235e3     7.9433e22     2e-2      -3   153e3   ;   % Wet Anorthite60 (Rybacki and Dresen, 2004)
           11   2.5119e15      1e11     2800      30      20       10e6         1e6      0.4     2e-7    2e-5   45e-13      1e3       273      1.6      3.0  345e3     1.9953e22     2e-2      -3   159e3   ;   % Wet Anorthite (Rybacki et al., 2006; siehe visc calc)
           12   2.5119e15      1e11     2800      30      20       10e6         1e6      0.4     2e-7    2e-5   45e-13      1e3       273      1.6      3.0  345e3     1.9953e22     2e-2      -3   159e3   ;   % Wet Anorthite (Rybacki et al., 2006; siehe visc calc)
           13   4.6627e15      1e11     3000      30      20       1e6          1e6      1.0     2e-7  2.5e-5   45e-13      1e3       273      2.0      4.1  723e3     8.2645e18     1e-3      -3   436e3   ;   % Oceanic Crust; An50Di50 (Diamov and Dresen, 2005)
           14   4.6627e15      1e11     3000      30      20       10e6         1e6      0.0     2e-7  2.5e-5   45e-13      1e3       273      2.0      4.1  723e3     8.2645e18     1e-2      -3   436e3   ;   % Oceanic Crust; An50Di50 (Diamov and Dresen, 2005)
           15   9.0909e15      1e11     3300      30      20       10e6         1e6      0.0     2e-8  2.5e-5   45e-13      1e3       273      1.6      3.5  480e3     6.6667e14     1e-3   -3.00   335e3   ;   % Upper Mantle; Wet olivine (Hirth and Kohlstedt, 2003)
           16   9.0909e15      1e11     3300      30      20       10e6         1e6      0.0     2e-8  2.5e-5   45e-13      1e3       273      1.6      3.5  480e3     6.6667e14     1e-3   -3.00   335e3   ;   % Upper Mantle; Wet olivine (Hirth and Kohlstedt, 2003)
           17   9.0909e15      1e11     3300      30      20       10e6         1e6      0.0     2e-8  2.5e-5   45e-13      1e3       273      1.6      3.5  480e3     6.6667e14     1e-3   -3.00   335e3   ;   % Incoming Mantle; Wet olivine (Hirth and Kohlstedt, 2003)
           18   9.0909e15      1e11     3300      30      20       10e6         1e6      0.0     2e-8  2.5e-5   45e-13      1e3       273      1.6      3.5  480e3     6.6667e14     1e-3   -3.00   335e3   ;   % Incoming Mantle; Wet olivine (Hirth and Kohlstedt, 2003)
%            16   9.0909e15      1e11     3300      30      30       25e6        10e6      0.0     2e-8    2e-5   45e-13      1e3       273      1.6      3.5  530e3     6.6667e14     1e-2   -3.00   375e3   ;   % Upper Mantle; Dry Olivine (Hirth and Kohlstedt, 2003)
%            17   9.0909e15      1e11     3300      30      30       25e6        10e6      0.0     2e-8    2e-5   45e-13      1e3       273      1.6      3.5  530e3     6.6667e14     1e-2   -3.00   375e3   ;   % Incoming Mantle; Dry Olivine (Hirth and Kohlstedt, 2003)
%            15   5.2830e49      1e11     3300      30      30       25e6        10e6      0.0     2e-8    2e-5   45e-13      1e3       273      1.6      8.2  682e3     1.2559e16    50e-6   -3.00   484e3   ;   % Upper Mantle; Olivine (Faul and Jackson, 2007; Faul et al., 2011)
%            16   5.2830e49      1e11     3300      30      30       25e6        10e6      0.0     2e-8    2e-5   45e-13      1e3       273      1.6      8.2  682e3     1.2559e16    50e-6   -3.00   484e3   ;   % Upper Mantle; Olivine (Faul and Jackson, 2007; Faul et al., 2011)
%            17   5.2830e49      1e11     3300      30      30       25e6        10e6      0.0     2e-8    2e-5   45e-13      1e3       273      1.6      8.2  682e3     1.2559e16    50e-6   -3.00   484e3   ;   % Incoming Mantle; Olivine (Faul and Jackson, 2007; Faul et al., 2011)
%            15   3.5481e23      1e11     3300      30      30       25e6        10e6      0.0     2e-8    2e-5   45e-13      1e3       273      1.6     4.94  610e3     4.2658e18        1   -2.98   261e3   ;   % Upper Mantle; Olivine (Korenang and Karato, 2008)
%            16   3.5481e23      1e11     3300      30      30       25e6        10e6      0.0     2e-8    2e-5   45e-13      1e3       273      1.6     4.94  610e3     4.2658e18        1   -2.98   261e3   ;   % Upper Mantle; Olivine (Korenang and Karato, 2008)
%            17   3.5481e23      1e11     3300      30      30       25e6        10e6      0.0     2e-8    2e-5   45e-13      1e3       273      1.6     4.94  610e3     4.2658e18        1   -2.98   261e3   ;   % Incoming Mantle; Olivine (Korenang and Karato, 2008)
%            18    4.88e133      1e11     2800      30      30        1e7         1e7      0.0     2e-6    2e-5   45e-13      1e3       273      2.5     18.0   51e3             0        0       0   220e3   ;   % Biotite (Kronenberg et al., 1990)
           19   9.0909e15      1e11     3300      15      15        1e6         1e6      1.0     2e-8  2.5e-5   45e-13      1e3       273      1.6      3.5  480e3     6.6667e14     1e-3   -3.00   335e3   ;   % Upper Mantle; Wet olivine (Hirth and Kohlstedt, 2003)
           20    4.88e133      1e11     2800      30      30        1e7         1e7      0.0     2e-6    2e-5   45e-13      1e3       273      2.5     18.0   51e3             0        0       0   220e3   ];   % Biotite (Kronenberg et al., 1990)

                 
Phase   =   ROCKS(:,1)';
eta     =   ROCKS(:,2)';
mu      =   ROCKS(:,3)';
rho     =   ROCKS(:,4)';
phi     =   ROCKS(:,5)';
phi_w   =   ROCKS(:,6)';
C       =   ROCKS(:,7)';
C_w     =   ROCKS(:,8)';
lambda  =   ROCKS(:,9)';
if Temperature == 1
    hr  =   ROCKS(:,10)';
    ta  =   ROCKS(:,11)';
    tb  =   ROCKS(:,12)';
    cp  =   ROCKS(:,13)';
    T   =   ROCKS(:,14)';
    kk  =   ROCKS(:,15)';
end
if Powerlaw == 1
    n   =   ROCKS(:,16);
    Q   =   ROCKS(:,17);
end
if Diffusion_creep == 1
    dCr =   ROCKS(:,18);
    Gr  =   ROCKS(:,19);    
    m   =   ROCKS(:,20);
    Qd  =   ROCKS(:,21);
end
    
OH_const = 50e6;
Lith = 120e3;
Subd_vel = 0.01;
Plate_age = 75e6;

% ====================== variables to chose ===============================

lambda_linear   = 3;        % 1: initial lambda; 2: linear increase of lambda of 0.1 per 2000 m; 3: surface lambda is 0.4, bottom of sed lambda is 0.9, increases linearly between
hydr_circ_thick = 1e3;      % thickness of uppermost sediment pile prone to hydrostatic fluid pressure due to free fluid circulation
lambda_increase = 0.125;    % increase in lambda perkilometer in case of "lambda_linear = 2"
lambda_bottom   = 0.95;      % lambda at sediment base
lambda_dry_OC   = 0.4;      % lambda in dry oceanic crust
lambda_hyd_OC   = 0.99;     % lambda in hydrated oceanic crust
SediThick       = 1e3;      % Sediment thickness in meters

% ====================== variables to chose ===============================


% Geometry      Phase       x-min       x-max       y-min       y-max
SETUP   =   [   1           0           Lx          0           Ly      ;
                4           0           Lx       10e3         15e3      ;
                3           0           Lx     14.7e3         15e3      ];
      
% Plastic weakening thresholds
w_1 = 0.1;
w_2 = 1.0;

% Parameters
gravity_y       =   9.81;
gravity_x       =   0.0;
SecYear         =   3600*24*365.25;
p_init          =   1e3;

% Surface process
SedimentationStyle  =   2;          % 1 = sed and ero || 2 = only sed || 3 = only ero || 4 = linear below sealevel || 5 = sed and ero independent of waterlevel
SurfaceCoeff        =  [1e-6];   % Coefficient for [SEDIMENTATION EROSION] FOR DIFFERENT SURFACE PROCESSES, OTHERWISE ONLY ONE COEFFICIENT
SurfaceBC_left      =   -1;       % For free slip: -1
SurfaceBC_right     =    10e3;       % For free slip: -1
Surface_dt          =   1;         % makes surface process every xx timestep
time_old            =   0;
SediMarker          =  [ 7 8 ];       % Marker type for sedimentation
ErosMarker          =   1;          % Marker type for erosion
SedChange           =   1e6;        % Marker change for sedimentation in years
WaterLevel          =   0;          % If WaterLevel == 0 it is switched off !!
SediRate            =   0;      % Linear sedimentation rate in m/yr
surface_nodes       =   5;          % Times x nodes along surface
surface_init        =   10000;
surface_x           =   [0:(Lx/(nx-1))/surface_nodes:Lx];
surface_y           =   [surface_init*ones(1,length(surface_x))];
surface_smoother    =   3;


% Cutoff viscosities
eta_max = 1e24;
eta_min = 1e17;
% eta_bingham = 5e18;

if Temperature==1
    % Temperature boundary conditions
    A_thick     =   10e3; % Thickness of sticky-air
    C_thick     =   30e3; % Thickness of Crust
    L_thick     =   Lith; % Thickness of Lithosphere
    Moho_temp   =   500;  % Temperature at Moho in ?C
    L_A_temp    =   1270+Lith/2000; % Temperature at Lithosphere/Asthenosphere boundary in ?C
    M_grad      =   0.5;  % Temperature gradient in Mantle (background)
    Pert_beg    =   450e3;
    Pert_end    =   550e3;
    Pert_add    =   0e3; % km difference for Lithosphere/Asthenosphere depth
    top_T       =   273;    % temperature in case of shear model (in Kelvin)
    Ext_T       =   0.995;
    if Ext_T == 1
        bottom_T    =   L_A_temp+273+M_grad*(Ly-A_thick-L_thick)/1e3;
    elseif Ext_T == 2
        bottom_T    =   top_T;
    else
        bottom_T    =   ((L_A_temp+273)+M_grad*(Ly-L_thick-A_thick)/1000+(2000)/1000*M_grad/2) - Ext_T*((L_A_temp+273)+M_grad*(Ly-L_thick-A_thick)/1000-(2000)/1000*M_grad/2);
    end
    shear_heat  =   0.99;
    ra_heat     =   1;
    adiab_heat  =   1;
    TPdep_dens  =   1;
end

% Velocity boundary conditions

%   L T T T T T T T T T T T RR      T T T T T T T T T T T T T   
%   L                       RR      L                       R
%   L      x-velocity       RR      L      y-velocity       R
%   L                       RR      L                       R
%   L                       RR      B B B B B B B B B B B B B
%   L B B B B B B B B B B B RR      B B B B B B B B B B B B B


BC_left_init     =  [zeros(1,ny-2) -0.0025 -0.0075 -0.0125];  % no need if mirrored
BC_right_init    =  -0.01;
BC_top_init      =  -0.01*14.9/100;
BC_bottom_init   = -0.01;
% BC_bottom_init   =(-0.00*135/1000)-0.99*(-0.00*132.989/1000);

bound            = {'freeslip','noslip','velocity','external','mixed','mirror'};
top_BC           = bound(1);
bottom_BC        = bound(3);
left_BC          = bound(2);
right_BC         = bound(1);
left_mix         = -ones(1,ny-2);
right_mix        = -ones(1,ny-2);
left_mix(20:110)  = 1;
right_mix(20:110) = 1;
% top_mix          = ones(1,nx-2);
bottom_mix(1,:)    = ones(1,nx-2);
bottom_mix(2,:)    = -0.005:0.01/(nx-3):0.005;
Ext_BC           = 0.99;

% Inversion
inversion = [100e6 120e6 135.5e6 136e6];     % Time in Ma for velocity inversion

% Incoming marker type
mark_top    = 1;
mark_bottom = 17;

%=============================================
% Initialize matrices and coordinates
initialize_beg;
%=============================================

% Initial marker pattern
pattern_marker  =   [ 4 ];     % Types of markers with pattern (check marker phases !!!!)
pattern_type    =   {'horizontal','vertical','kaki'};   % Stripes or kaki
pattern_type    =   pattern_type(3);
pattern_xdim    =   0.5e3;    % thickness of vertical stripes
pattern_ydim    =   0.5e3;    % thickness of horizontal stripes

%=============================================
% Initialize phase and temperature distribution of marker 
marker_distribution;
%=============================================

% Initial time and iteration setup
time        = 0;
t_beg       = 1;
dt_value    = 0.25;          % threshold value for maximal timestep, 0.1 = 10% movement of dx or dy
Ddt         = 1000*SecYear;  % Initial timestep to pre-stress
short_dt    = 1000*SecYear;  % Short initial timestep
dt_max      = 1000*SecYear;  % Timestep
n_short_dt  = 1;            % Number of short initial timesteps
miniter     = 1;             % minimal number of iterations
maxiter     = 1;             % maximal number of iterations
error       = 1;             % 1 = Error 1 (average nodal velocity change); 2 = Error 2 (largest nodal velocity change)
vel_res     = 1e-14;         % sum(velocity) change fot iter break
strain_rate_smoother = 0.87;

% Load breakpoint file if necessary
if exist(char('Breakpoint.mat'))
    load Breakpoint.mat
    t_beg    = timestep+1;
end

%=========================== START TIME LOOP ==========================
%======================================================================

for timestep = t_beg:3000
    tic

    dt=dt_max;
    if timestep <= n_short_dt
        dt = short_dt;  
    end
    
    velocity_inversion;
    
    for niter = 1:maxiter

        fprintf('Timestep: %d\n', timestep);
        fprintf('Iteration: %d\n', niter);
        
        %=============================================
        % Reload old values and initialize matrices
        initialize_iter;
        %=============================================
        
        %=============================================
        % Calculate power-law and brittle viscosities
        viscosity_calculation;
        %=============================================
        
        %=============================================
        % Fill nodal values from marker information
        marker_to_nodes; % dispatched loops, better for desktop
        toc, fprintf('for updating nodal values!\n'); 
        %=============================================   
        
        if Temperature == 1            
            % Boundary conditions to interpolated T
            % Upper BC
            Temp((ny+1)+1:(ny+1):(nx-1)*(ny+1)+1)   =   top_T*2 - Temp((ny+1)+2:(ny+1):(nx-1)*(ny+1)+2);
            % Lower BC
            if Ext_T == 1 || Ext_T == 2
                Temp(2*(ny+1):(ny+1):(nx)*(ny+1))       =   bottom_T'*2 - Temp(2*(ny)+1:(ny+1):(nx)*(ny+1)-1);
            else
                Temp(2*(ny+1):(ny+1):(nx)*(ny+1))       =   bottom_T + Ext_T*Temp(2*(ny)+1:(ny+1):(nx)*(ny+1)-1);
            end
            % Left BC
            Temp(1:ny+1)                            =   Temp(ny+2:2*(ny+1));
            % Right BC
            Temp(nx*(ny+1)+1:(nx+1)*(ny+1))         =   Temp((nx-1)*(ny+1)+1:(nx)*(ny+1));
        end
        
        %=====================================================
        
        %=============================================
        % Fill and solve matrix for stokes
        stokes_direct_solver;
        toc, fprintf('for solving the matrix!\n')
        %=============================================
        
        Vx    =   S(1:(ny+1)*(nx+1));
        Vy    =   S(1+(ny+1)*(nx+1):2*(ny+1)*(nx+1));
        P     =   S(2*(ny+1)*(nx+1)+1:end).*kcont;

        P(1:ny+1)         =   P(ny+2:(ny+1)*2);
        P(end-ny:end)     =   P(end-(ny+1)-ny:end-(ny+1));
        P(1:ny+1:end-ny)  =   P(2:(ny+1):end-ny+1);
        P(ny+1:ny+1:end)  =   2.*P(ny:ny+1:end-1)-P(ny-1:ny+1:end-2);

        
        Vx2d    =   reshape(Vx,ny+1,nx+1);
        Vy2d    =   reshape(Vy,ny+1,nx+1);
        P2d     =   reshape(P,ny+1,nx+1);
                
        %======================================================================
        % redistributing S
        
        %   1--6--11--16
        %   |  |  |    |
        %   2--7--12--17
        %   |  |  |    |
        %   3--8--13--18
        %   |  |  |    |
        %   4--9--14--19
        %   |  |  |    |
        %   5-10--15--20
        
        % Vxs      =   (Vx(1:end-(ny+1)) .*dy_2d(2:end-(ny+1)+1) + Vx(2:end-(ny+1)+1).*dy_2d(1:end-(ny+1)))./(dy_2d(2:end-(ny+1)+1) + dy_2d(1:end-(ny+1))); 
        % Vys      =   (Vy(1:end-(ny+1)) .*dx_2d(ny+2:end) + Vy(ny+2:end).*dx_2d(1:end-(ny+1)))./(dx_2d(ny+2:end) + dx_2d(1:end-(ny+1)));
                
        % Vx2d_p(:,2:nx+1)   =   (Vx2d(:,1:end-1) + Vx2d(:,2:end))./2;
        % Vy2d_p(2:ny+1,:)   =   (Vy2d(1:end-1,:) + Vy2d(2:end,:))./2;
        % 
        % Vx2d_p(:,1)        =   -Vx2d_p(:,2);
        % Vx2d_p(:,end)      =   -Vx2d_p(:,end-1);
        % Vy2d_p(1,:)        =   -Vy2d_p(2,:);
        % 
        % Vxp    =   reshape(Vx2d_p,(ny+1)*(nx+1),1);
        % Vyp    =   reshape(Vy2d_p,(ny+1)*(nx+1),1);
        % 
        % Vxs((ny+1)*nx:(ny+1)*(nx+1))=0;
        % Vxs(ny+1:ny+1:end)=0;
        % Vys((ny+1)*nx:(ny+1)*(nx+1))=0;
        % Vys(ny+1:ny+1:end)=0;
        % 
        % define optimal timestep
        Ddt = dt;
        
        Vx_max  =   max(abs(Vx));
        Vy_max  =   max(abs(Vy));
        
        if dt_value*1/(Vx_max/(dx) + Vy_max/(dy)) < Ddt
            Ddt  =   dt_value*1/(Vx_max/((dx)) + Vy_max/(dy));
        end
                
        %=============================================
        % Calculating strain rates and stresses
        strain_rate_stress_calculation;
        toc, fprintf('for strain/stress calculation and interpolation to markers!\n');
             fprintf('========= Ddt = %d years =============================\n',Ddt/SecYear);
        %=============================================
        
        fprintf('========= Min Pressure = %d MPa ===========================\n',min(P)/1e6);
        fprintf('========= Max Pressure = %d MPa ===========================\n',max(P)/1e6);
   
        %=============================================
        % Calculating strain rates and stresses
        vel_res_calculation;
        if timestep>1
            fprintf('========= VELOCITY ERROR %d = %d =========\n\n',error,iterations(error,niter+maxiter*(timestep-1)));
            fprintf('========= Pressure change on node = %d =========\n\n',dP_node(timestep));
        end
        %=============================================

        % Plotting if running on desktop
        if brutus==0 && timestep>1
            E2nd_s_2d   =   reshape(E2nd_s,ny+1,nx+1);
            eta_s_2d    =   reshape(eta_s,ny+1,nx+1);
            strain_2d   =   reshape(strain,ny+1,nx+1);
            strainv_2d   =   reshape(strainv,ny+1,nx+1);
            lambda_s_2d   =   reshape(lambda_s,ny+1,nx+1);
            
            figure(1), clf
            colormap jet
            
            subplot(231)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),log10(E2nd_s_2d(1:end-1,1:end-1)))%, caxis([-14 -10])
            shading interp
            colorbar
            axis image, axis ij
            title('E2nd [1/s]')
            
            subplot(232)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),log10(eta_s_2d(1:end-1,1:end-1))), caxis(log10([eta_min eta_max]))
            shading interp
            colorbar
            axis image, axis ij
            title(['\eta_{s} [Pa.s] ',    num2str(Ddt/SecYear)])
            drawnow

            subplot(233)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),(strain_2d(1:end-1,1:end-1)))
            shading interp
            colorbar
            axis image, axis ij
            title(['Strain [-] ',    num2str(Ddt/SecYear)])
            drawnow
            
            subplot(234)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),(T2nd_s_2d(1:end-1,1:end-1)))
            shading interp
            colorbar
            axis image, axis ij
            title(['\tau_{s} [Pa] ',    num2str(Ddt/SecYear)])
            drawnow

            subplot(235)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),Vx2d_s*SecYear)
            shading interp
            colorbar
            axis image, axis ij
            title(['V_{x} [m/yr] ',    num2str(Ddt/SecYear)])
            drawnow

            subplot(236)
            pcolor(x_Vx(1:end-1),y_Vy(1:end-1),Vy2d_s*SecYear)
            shading interp
            colorbar
            axis image, axis ij
            title(['V_{y} [m/yr] ',    num2str(Ddt/SecYear)])
            drawnow

            
            figure(444)
            subplot(211), semilogy((niter-1)/maxiter+(timestep-1),iterations(1,niter+maxiter*(timestep-1)),'ko'), hold on
            title('Error 1: with sticky-air')
            subplot(212), semilogy((niter-1)/maxiter+(timestep-1),iterations(2,niter+maxiter*(timestep-1)),'ko'), hold on
            title('Error 2: with sticky-air')
        end
        
        % Exiting the iteration loop
        if (niter>=miniter && iterations(error,niter+maxiter*(timestep-1))<vel_res) || timestep==1
            break;
        end
    end

    %==========================================================================
    dt = Ddt;
    
    % Calculating accumulated plastic/viscous strain and grain size
    Strain_GrainSize_GSE;
    %==========================================================================

    %=============================================
    % Calculating temperature
    subgrid_diffusion_stresses;
    stress_rotation;
    % Average deviatoric stress
%     Sigma_d(1,timestep)   =   sum((Txxm.^2 + Txym.^2).^0.5)./length(T2ndm);
%     Visc_d(1,timestep)   =   sum(eta_effm)./length(eta_effm);
    %=============================================
    
    %=============================================
    % Calculating temperature
    if Temperature == 1
        temperature_direct_solver;
        ind=Tm<273;
        Tm(ind)=273;
        toc, fprintf('for temperature calculation!\n');
    end
    %=============================================


    %=============================================
    % Move markers
    move_marker;
    move_surface;
    if Beam_function == 1
        beam_function;
    end

    toc, fprintf('for moving makers!\n');
    %=============================================

    time = time + dt;

    %======================================================================
    % Visualization begin
    %======================================================================


    if brutus == 0  figure(3), clf
        
        for i=1:20
            Phase1  =   find(Im == i);
            marksize = 5;
            hold on
            plot(xm(Phase1),ym(Phase1),'.','Color',Colormaps(i,:),'MarkerSize',marksize)
        end
      
        
%         plot(surface_x,surface_y,'r'), hold on
        
        axis image, axis ij
        title(['Time: ',num2str((time)/SecYear/1e6),' Ma'])
%         plot(x2d_b,y2d_b,'k',x2d_b',y2d_b','k')
%         quiver(x2d_b,y2d_b,Vx2db,Vy2db,'r','LineWidth',1.5)
        E2nd_s_2d   =   reshape(E2nd_s,ny+1,nx+1);
        eta_s_2d    =   reshape(eta_s,ny+1,nx+1);
        T2nd_s_2d   =   reshape(T2nd_s,ny+1,nx+1);
                
        hold off
                
%         figure(8)
%         plot(surface_x,surface_y,'b'), hold off
%         axis ij
        
        figure(4), clf
        colormap jet
        
        subplot(311)
        pcolor(x_Vx(1:end-1),y_Vy(1:end-1),log10(E2nd_s_2d(1:end-1,1:end-1)))%, caxis([-14 -10])
        shading interp
        colorbar
        axis image, axis ij
        title('E2nd ')
        
        subplot(312)
        pcolor(x_Vx(1:end-1),y_Vy(1:end-1),log10(eta_s_2d(1:end-1,1:end-1))), caxis(log10([eta_min eta_max]))
        shading interp
        colorbar
        axis image, axis ij
        title(['\eta_{s}   ',    num2str(dt/SecYear)])
        
        subplot(313)
        pcolor(x_Vx(1:end-1),y_Vy(1:end-1),(T2nd_s_2d(1:end-1,1:end-1)))
        shading interp
        colorbar
        axis image, axis ij
        title(['T2nd'])
        drawnow
    end

    %==========================================================================
    % Visualization end
    %==========================================================================

    %=============================================
    % Delete markers out of grid
    outgrid_marker;
    %=============================================

    %=============================================
    % New markers from sides
    incoming_marker;
    toc, fprintf('for outgoing/incoming markers!\n');
    %=============================================
    
    %=============================================
    % Surface process
    if  Surface_diffusion == 1
        surface_calculation;
        toc, fprintf('for surface process!\n');
    end
    %=============================================
    work_ratem = Txxm*2.*Exxn + Txym.*2.*Exyn;
    workm = workm + work_ratem*dt;
    

    fprintf('dt   = %d years\n', dt/SecYear)
    fprintf('Time = %d years\n\n', time/SecYear)

    toc, fprintf('for complete timestep\n===============\nMarknum: %d\n===============\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n', marknum);    

    %======================== SAVE OUTPUT FILE ============================
    if mod(timestep,create_output)==0 || timestep==2 || timestep==3
        
        fname = ['Shearing_',num2str(1e6+timestep),'.mat'];
        cd OUTPUT
        save([char(fname)]   ,'Vx2d','Vy2d','E2nd_s','T2nd_s','xm','Txx','Txy',...
            'ym','Im','eta_s','P2d','Lx','Ly','strain','strainv','lambda_s','rho_s','surface_x','surface_y',...
            'nx','ny','time','dt','SecYear','x','y','iterations','Pm','workm','work')
        cd ..
    end

    % Save breakpoint file
    if mod(timestep,create_breakpoint)==0
        fname = ['Breakpoint1.mat'];
        save(char(fname))
%         break
        movefile('Breakpoint1.mat','Breakpoint.mat');
    end

end
