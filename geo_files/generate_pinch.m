function [diameter] = generate_pinch(N,aspect,pamp,n,plotit)
%
% Racetrack with periscribed pinch peristalsis
%
% Will return the diameter of the tube for the following inputs: 
%
%   N = number of Cartesian grid meshwidths at the finest level of the AMR grid 
%   aspect = the length to diameter aspect ratio of the tube (the length of
%       the tube, Let, is going to stay the same.
%   pamp = the compression ratio or percent occlusion of the tube during a
%       beat
%   plotit = logical indicating wish to plot geometry variables (1) or not
%   (0 or anything else)
%
% Example use: [dia]=generate_pinch(512,10,0.95,1,1) where 512 is the grid
% size, 10 is the aspect ratio of the tube, 0.95 is the compression ratio,
% and 1 is the run number associated with these parameter values.
%

% Parameters for the IBAMR input2d setup
L = 1;          % Length of computational domain (m)
%N = 512;       % Assigned by user.
dx = L/N;       % Cartesian mesh width (m)
ds = L/(2*N);   % space between boundary points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for racetrack

Let = 0.4;              % Length of elastic section of tube (m)
Nend = 10;              % Number of rigid points on each end of the elastic section
Lt = Let + 2*Nend*ds;   % Length of straight section with rigid points on each end. 
%aspect = 4;            % Assigned by user.

diameter = Let/aspect         % diameter of the tube (m)
R2 = 0.1;               % radius of the inner wall
R1 = R2 + diameter;     % radius of the outer wall

Nstraight = 2*ceil(Lt/ds);  % number of points along each straight section
Ncurve = 2*ceil(pi*R1/ds);  % number of points along each curved section
Nrace = Nstraight+2*Ncurve; % number of points making up the racetrack. 
dtheta = pi/(Ncurve/2);     % angle increment for drawing curved edges.

mesh_name = 'heart_';       % structure name

centery = 0;            % y-position of center of curved sections
centerx1 = -0.5*Lt;     % x-position of center of left curved section
centerx2 = 0.5*Lt;      % x-position of center of right curved secton

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for perscribed peristaslsis state changes
sigma = 0.005;          % Determines pointiness of pinch
mu = 0.19;              % Determines how far from the x-center the pinch starts
%pamp = 0.8;            % Assigned by user.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for markers
Nmarkersx = 11;                     %number of columns of markers
Nmarkersy = 11;                     %number of markers in each column
Nmarkers=Nmarkersx*Nmarkersy;       %total number of markers
dmx = Let/(Nmarkersx-1);            %space between markers in x-direction
dmy = diameter/(Nmarkersy-1);       %space between markers in y-direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% material parameters
kappa_spring = 30.0;               % spring constant (Newton)
kappa_beam = 0.3;                 % beam stiffness constant (Newton m^2) %2.5e-2 works for Wo>=5
kappa_target = kappa_spring;        % target point penalty spring constant (Newton)
Fmag = 4.0e0;                % this is my best guess at a reasonable applied force %4.0e0 works for Wo>=5
phase = 0;                      %initial phase of the oscillating force, where F=Fmag*phase and phase = (1+sin(2*pi*f*t-pi/2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if plotit == 1 
% Initialize plotting function
    figure(1) 
    hold on
    ylim([-L/2 L/2])
    xlim([-L/2 L/2])
else
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out elastic section of tube

% Allocate space for variables
ytop_elastic = zeros(1,ceil(Nstraight/2));
xtop_elastic = zeros(1,ceil(Nstraight/2));
ybot_elastic = zeros(1,ceil(Nstraight/2));
xbot_elastic = zeros(1,ceil(Nstraight/2));

% Vertex information
vertex_fid = fopen([mesh_name 'tube_' num2str(n) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', Nstraight);

% Top section, elastic tube
for i = 1:ceil(Nstraight/2)
    ytop_elastic(1,i) = centery-R2;
    xtop_elastic(1,i) = -Lt/2+(i-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop_elastic(1,i), ytop_elastic(1,i));
end

% Bottom section, elastic tube
for  i = 1:ceil(Nstraight/2)
    ybot_elastic(1,i) = centery-R1;
    xbot_elastic(1,i) = -Lt/2+(i-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot_elastic(1,i), ybot_elastic(1,i));
end

fclose(vertex_fid);

if plotit==1
% Plots elastic tube vertices
    plot(xtop_elastic,ytop_elastic,'r.')
    plot(xbot_elastic,ybot_elastic,'r.')
else 
    
end
% Write out the spring information for the elastic section

spring_fid = fopen([mesh_name 'tube_' num2str(n) '.spring'], 'w');
fprintf(spring_fid, '%d\n', Nstraight-2);

%elastic part of tube
for i = 0:ceil(Nstraight/2)-2,
    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', i, i+1, kappa_spring*ds/(ds^2), ds);
end

for i = ceil(Nstraight/2):Nstraight-2,
    fprintf(spring_fid, '%d %d %1.16e %1.16e\n', i, i+1, kappa_spring*ds/(ds^2), ds);
end

fclose(spring_fid);

% Write out the beam information for the elastic section

beam_fid = fopen([mesh_name 'tube_' num2str(n) '.beam'], 'w');
fprintf(beam_fid, '%d\n', Nstraight-4);

%elastic part of tube
for i = 0:ceil(Nstraight/2)-3,
    fprintf(beam_fid, '%d %d %d %1.16e\n', i, i+1, i+2, kappa_beam*ds/(ds^4));
end

for i = ceil(Nstraight/2):Nstraight-3,
    fprintf(beam_fid, '%d %d %d %1.16e\n', i, i+1, i+2, kappa_beam*ds/(ds^4));
end
fclose(beam_fid);

% Write out the target point information for the ends of the elastic tube
target_fid = fopen([mesh_name 'tube_' num2str(n) '.target'], 'w');

fprintf(target_fid, '%d\n', 4*Nend);

for i = 0:Nend-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

for i = ceil(Nstraight/2)-Nend:ceil(Nstraight/2)-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

for i = ceil(Nstraight/2):ceil(Nstraight/2)+Nend-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

for i = Nstraight-Nend:Nstraight-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

fclose(target_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make markers as vertices with no material properties

% Allocates space
x_mark = zeros(Nmarkersx,Nmarkersy);
y_mark = zeros(Nmarkersx,Nmarkersy);

vertex_fid = fopen(['markers_' num2str(n) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', Nmarkers);

% Markers
for i=0:Nmarkersx-1,
    for j=0:Nmarkersy-1,
        y_mark(i+1,j+1) = centery-R2-j*dmy;
        x_mark(i+1,j+1) = -Let/2+i*dmx;
    fprintf(vertex_fid, '%1.16e %1.16e\n', x_mark(i+1,j+1), y_mark(i+1,j+1));
    end
end

if plotit==1
% Plot markers
    plot(x_mark,y_mark,'b.')
else
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% race track part

% Allocate Space
x_race = zeros(1,Nrace);
y_race = zeros(1,Nrace);

% Write out the vertex information
vertex_fid = fopen([mesh_name 'race_' num2str(n) '.vertex'], 'w');
fprintf(vertex_fid, '%d\n', Nrace);

%right curved part of racetrack
for i=1:ceil(Ncurve/2),
    theta = (i-1)*dtheta-pi/2;
    y_race(1,i) = centery+R2*sin(theta);
    x_race(1,i) = Lt/2+R2*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', x_race(1,i), y_race(1,i));
end

for i=ceil(Ncurve/2)+1:Ncurve,
    theta=(i-Ncurve/2-1)*dtheta-pi/2;
    y_race(1,i) = centery+R1*sin(theta);
    x_race(1,i) = Lt/2+R1*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', x_race(1,i), y_race(1,i));
end

%straight section on the top
for i = Ncurve+1:Ncurve+ceil(Nstraight/2),
    y_race(1,i) = centery+R2;
    x_race(1,i) = centerx2-(i-Ncurve-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', x_race(1,i), y_race(1,i));
end

for i = Ncurve+ceil(Nstraight/2)+1:Ncurve+Nstraight,
    y_race(1,i) = centery+R1;
    x_race(1,i) = centerx2-(i-Ncurve-ceil(Nstraight/2)-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', x_race(1,i), y_race(1,i));
end

%left curved part of racetrack
for i = Ncurve+Nstraight+1:Ncurve+Nstraight+ceil(Ncurve/2),
    theta = pi/2+(i-Ncurve-Nstraight-1)*dtheta;
    y_race(1,i) = centery+R2*sin(theta);
    x_race(1,i) = centerx1+R2*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', x_race(1,i), y_race(1,i));
end

for i = Ncurve+Nstraight+ceil(Ncurve/2)+1:2*Ncurve+Nstraight,
    theta = pi/2+(i-Ncurve-Nstraight-ceil(Ncurve/2)-1)*dtheta;
    y_race(1,i) = centery+R1*sin(theta);
    x_race(1,i) = centerx1+R1*cos(theta);
    fprintf(vertex_fid, '%1.16e %1.16e\n', x_race(1,i), y_race(1,i));
end
fclose(vertex_fid);

if plotit==1
% Plots the racetrack vertices
    plot(x_race,y_race,'k.')
else
    
end

% Write out the target point information for the racetrack
target_fid = fopen([mesh_name 'race_' num2str(n) '.target'], 'w');

fprintf(target_fid, '%d\n', Nrace);

for i = 0:Nrace-1,
    fprintf(target_fid, '%d %1.16e\n', i, kappa_target*ds/(ds^2));
end

fclose(target_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out state change files

% Allocate space
ytop_elastic_s1 = zeros(1,ceil(Nstraight/2));
xtop_elastic_s1 = zeros(1,ceil(Nstraight/2));
ybot_elastic_s1 = zeros(1,ceil(Nstraight/2));
xbot_elastic_s1 = zeros(1,ceil(Nstraight/2));
ytop_elastic_s2 = zeros(1,ceil(Nstraight/2));
xtop_elastic_s2 = zeros(1,ceil(Nstraight/2));
ybot_elastic_s2 = zeros(1,ceil(Nstraight/2));
xbot_elastic_s2 = zeros(1,ceil(Nstraight/2));
ytop_elastic_s3 = zeros(1,ceil(Nstraight/2));
xtop_elastic_s3 = zeros(1,ceil(Nstraight/2));
ybot_elastic_s3 = zeros(1,ceil(Nstraight/2));
xbot_elastic_s3 = zeros(1,ceil(Nstraight/2));

% Vertex information State 1 (straight lines)
vertex_fid = fopen([mesh_name 'tube_state1_' num2str(n) '.txt'], 'w');
%fprintf(vertex_fid, '%d\n', Nstraight);

% Top section, elastic tube
for i = 1:ceil(Nstraight/2)
    ytop_elastic_s1(1,i) = centery-R2;
    xtop_elastic_s1(1,i) = -Lt/2+(i-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop_elastic_s1(1,i), ytop_elastic_s1(1,i));
end

% Bottom section, elastic tube
for  i = 1:ceil(Nstraight/2)
    ybot_elastic_s1(1,i) = centery-R1;
    xbot_elastic_s1(1,i) = -Lt/2+(i-1)*ds;
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot_elastic_s1(1,i), ybot_elastic_s1(1,i));
end

fclose(vertex_fid);

% Vertex information State 2 (pinch at left end)
vertex_fid = fopen([mesh_name 'tube_state2_' num2str(n) '.txt'], 'w');
%fprintf(vertex_fid, '%d\n', Nstraight);

% Top section, elastic tube
for i = 1:ceil(Nstraight/2)
    xtop_elastic_s2(1,i) = -Lt/2+(i-1)*ds;
    ytop_elastic_s2(1,i) = centery-R2-(diameter*pamp/2)*exp(-0.5*((xtop_elastic_s2(1,i)+mu)/sigma).^2);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xtop_elastic_s2(1,i), ytop_elastic_s2(1,i));
end

% Bottom section, elastic tube
for  i = 1:ceil(Nstraight/2)
    xbot_elastic_s2(1,i) = -Lt/2+(i-1)*ds;
    ybot_elastic_s2(1,i) = centery-R1+(diameter*pamp/2)*exp(-0.5*((xbot_elastic_s2(1,i)+mu)/sigma).^2);
    fprintf(vertex_fid, '%1.16e %1.16e\n', xbot_elastic_s2(1,i), ybot_elastic_s2(1,i));
end

fclose(vertex_fid);

if plotit ==1
    plot(xtop_elastic_s2,ytop_elastic_s2,'m.')
    plot(xbot_elastic_s2,ybot_elastic_s2,'m.')
else
    
end

% Vertex information State 3 (pinch at right end)
vertex_fid = fopen([mesh_name 'tube_state3_' num2str(n) '.txt'], 'w');
%fprintf(vertex_fid, '%d\n', Nstraight);

% Top section, elastic tube
for i = 1:ceil(Nstraight/2)
    xtop_elastic_s3(1,i) = -Lt/2+(i-1)*ds;
    ytop_elastic_s3(1,i) = centery-R2-(diameter*pamp/2)*exp(-0.5*((xtop_elastic_s2(1,i)-mu)/sigma).^2);
    fprintf(vertex_fid, '%1.16e \t %1.16e\n', xtop_elastic_s3(1,i), ytop_elastic_s3(1,i));
end

% Bottom section, elastic tube
for  i = 1:ceil(Nstraight/2)
    xbot_elastic_s3(1,i) = -Lt/2+(i-1)*ds;
    ybot_elastic_s3(1,i) = centery-R1+(diameter*pamp/2)*exp(-0.5*((xbot_elastic_s2(1,i)-mu)/sigma).^2);
    fprintf(vertex_fid, '%1.16e \t %1.16e\n', xbot_elastic_s3(1,i), ybot_elastic_s3(1,i));
end

fclose(vertex_fid);

if plotit ==1
    plot(xtop_elastic_s2,ytop_elastic_s3,'b.')
    plot(xbot_elastic_s2,ybot_elastic_s3,'b.')
    
    hold off
else
    
end


