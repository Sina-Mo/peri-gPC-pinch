function please_Create_Geometry()


%
% 1. Function will create the ".vertex" and ".target" point files as
%    necessary inputs for an IBAMR simulation. 
% 2. Creates two different phase files (xP1.txt, yP1.txt, xP2.txt, yP2.txt)
%    and prints their coordinates to .txt files -> for reading into
%    update_target_point_positions.C
%


%
% Grid Paramater Information %
%
L = 2.0;                    %Length of Computational Domain [-L/2,L/2]
NFINEST = 1024;             %FINEST GRID RESOLUTION, i.e., L/1025 Eulerian spacing
ds = L/(2*NFINEST);         %Lagrangian Pt. Spacing


%
% Structure Parameter Information
%
r = 0.2;                    %Radii for circle geometry
geo_name = 'circle';        %Name of geometry for input files (.vertex, .target, etc)



%
% Calls a function that will create the circle (elliptical) geometry and
%         store as vectors, X and Y, for xPts and yPts respectively.
%
[X Y] = please_Create_Circle_Geometry(ds,r);



%
%Shift Geometry Accordingly (make two states) -> Circle centered at (-L/2,-L/2) and(L/2,L/2)
%      ->NOTE: For these state-type interpolations, you MUST make sure to have the same 
%              number of points in each state!
%
X1 = X - L/4;    Y1 = Y - L/4;  %Shifts one circle to left and down
X2 = X + L/4;    Y2 = Y + L/4;  %Shifts other circle to right and up


%
% PLOTS the Geometry you created above
%
plot(X1,Y1,'*'); hold on;
plot(X2,Y2,'r*'); hold on;
xlabel('x'); ylabel('y'); title('Structure Geometry:  Blue (Phase1)  &  Red (Phase2)'); axis([-L/2 L/2 -L/2 L/2]);


% 
% PRINT THE (X,Y)-coordinates into their own .txt file that will be read
%       in by update_target_point_position.C to interpolate between the states
% 
please_Print_Interpolation_State_Files(X1,Y1,X2,Y2)



%
% PRINT THE ACTUAL .VERTEX and .TARGET INFO TO FILES
%
target_force = 1e3;                                %Target Pt. "force"
please_Print_Vertex_File(X1,Y1,geo_name);          %Prints .vertex file (MAKE SURE IN 1ST STATE)    
please_Print_Target_File(X,geo_name,target_force); %Prints .target file





% Friendly Reminders! %
fprintf('\n\n\nThis function creates two input files: %s.vertex and %s.target\n',geo_name,geo_name);
fprintf('It also creates state geometry files in the form of .txt files\n\n\n'); 
fprintf('Don"t forget to copy .target & .vertex files to the other folder!\n\n\n');
fprintf('Remember to copy all the .txt files as well!\n\n\n');
fprintf('In update_target_point_positions, make sure "numPts=%d"\n\n\n',length(X1));
fprintf('Oh, and make sure the "L" value (L=%d) matches the L in input2d\n\n\n',L);
fprintf('Oh...one more thing, make sure your geometry names match those in input2d!\n');
fprintf('      (line 88 and 89 in input2d)\n\n\n'); 
fprintf('Eh, I guess check the update_target_point_positions.C file to make sure the\n');
fprintf('      interpolation makes sense\n\n\n');







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: creates CIRCLE (structure) geometry file 
% Returns: (X,Y)-pt. vectors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,Y] = please_Create_Circle_Geometry(ds,r)

% ds: Lagrangian spacing
% r: radii for circular geometry

dANG = ds/r;       % Incremental change in angle between points (used s = r*ANG)
ang = 0:dANG:2*pi; % Vector of angle values around circle

for i=1:length(ang)        % Loops over every angle in the angleVector
   X(i) = r*cos( ang(i) ); % Gives x-Values around Circle
   Y(i) = r*sin( ang(i) ); % Gives y-Values around Circle
end










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: prints the vertex pts to the .vertex file format
% CREATES: filename called "geo_name.vertex" (whatever geo_name is)
% RETURNS: null
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function please_Print_Vertex_File(X,Y,geo_name)            

% (X,Y): vectors containing each (x,y)-pt in the geometry
% geo_name: string for what you'd like file prefix to be called

Npts = length(X); %Gives # of pts. in geometry

vertex_fid = fopen([geo_name,'.vertex'], 'w');

% 1st line is the number of vertices in the file
fprintf(vertex_fid, '%d\n', Npts);

for i=1:Npts
        fprintf(vertex_fid, '%1.16e %1.16e\n',X(i), Y(i) );  %Prints in .vertex file format
end

fclose(vertex_fid);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: prints the target pts to the .target file format
% CREATES: filename called "geo_name.target" (whatever geo_name is)
% RETURNS: null
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function please_Print_Target_File(X,geo_name,target_force)


% (X,Y): vectors containing each (x,y)-pt in the geometry
% geo_name: string for what you'd like file prefix to be called
% target_force: target pt "force"

Npts = length(X); %Gives # of pts. in geometry

% Write out the target point information
target_fid = fopen([geo_name,'.target'], 'w');

fprintf(target_fid, '%d\n', Npts);

s = 0;
for n=1:Npts;
        fprintf(target_fid, '%d %1.16e\n', s , target_force);
        s = s+1;
end

fclose(target_fid);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: prints the X1, Y1, X2, and Y2 data to their own .txt file
% CREATES: xP1.txt, yP1.txt, xP2.txt, yP2.txt
% RETURNS: null
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function please_Print_Interpolation_State_Files(X1,Y1,X2,Y2)


fileID_xP1 = fopen('xP1.txt','w'); % Opens (creates) file called xP1.txt
fprintf(fileID_xP1,'%1.16e\n',X1); % Writes vector X1 to xP1.txt
fclose(fileID_xP1);                % Closes the file

fileID_yP1 = fopen('yP1.txt','w');
fprintf(fileID_yP1,'%1.16e\n',Y1);
fclose(fileID_yP1);

fileID_xP2 = fopen('xP2.txt','w');
fprintf(fileID_xP2,'%1.16e\n',X2);
fclose(fileID_xP2);

fileID_yP2 = fopen('yP2.txt','w');
fprintf(fileID_yP2,'%1.16e\n',Y2);
fclose(fileID_yP2);
