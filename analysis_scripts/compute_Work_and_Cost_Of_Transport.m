%-------------------------------------------------------------------------------------------------------------------%
%
% Author: Nicholas A. Battista
% Date: 08/08/2018
% Email:  nickabattista[at]gmail[.]com
% Current Institution: TCNJ
% Lab: Waldrop Lab, Miller Lab, Battista Lab
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Computes the Work & Cost of Transport for a specific simulation
%
% Inputs: None
%
% Outputs: 1. Cost of transport (COT_Sim)
%          2. Work done by tube for pumping (Work_Sim)
%
% Notes: 1. Set path to desired IBAMR hier_IB2d_data folder (line 50)
%        2. Starting index of desired data, e.g., F.01000 -> start = 1000 (line 43)
%        3. Ending index of desired data to analyze, e.g, F.28000 -> finish = 28000 (line 44)
%        4. The interval between successive data dumps, e.g., if F.01000 and next is F.02000, and next is F.03000 -> dump_int = 1000 (line 45)
%        5. Name of folder for new VTK tangential, magnitude of force data (line 51)
%        6. Set path to new VTK force data folder (line 56)
%        7. Flag to plot cost of transport and work over time (1=YES, 0=NO) (line 39)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [COT_Sim,Work_Sim] = compute_Work_and_Cost_Of_Transport(run,finish)

% Characteristic Length Scale
L = 0.1;            % Characteristic length scale (width of tube?) of system for COT  (line 129)

% TEMPORAL INFO FROM input2d %
dt = 1.0e-5;      % Time-step in simulation (from input2d)

% PLOT FLAG FOR COST OF TRANSPORT / WORK: yes/no
plot_Flag = 0;    % 1=yes, 0=no


% DATA ANALYSIS INFO %
start=20000;                              % 1ST interval # included in data analysis
%finish=250000;                            % LAST interval # included in data analysis 
dump_int = 5000;                         % Interval between successive data dumps for ANALYSIS
dump_Times = (start:dump_int:finish)*dt; % Time vector when data was printed for ANALYSIS analysis


% ANALYZE THE IBAMR DATA AND CONVERT TO .VTK

folder1 = ['hier_data_IB2d',num2str(run)];
path_name = ['/Volumes/HelmsDeep/IBAMR/peri-gPC/peri-pinch/hier_data/',folder1,'/']; % Path to IBAMR Data
folder_Name = 'forces_VTK_HELL_YEAH';                         % Name of folder to save .vtk data!
print_IBAMR_hier_data_to_VTK(path_name,start,dump_int,finish,dt,folder_Name);


% SET PATH TO NEW VTK FORCE DATA!
pathForce= ['/Volumes/HelmsDeep/IBAMR/peri-gPC/peri-pinch/hier_data/',folder1,'/',folder_Name,'/'];


% Initialize storage for force data 
fMagAvg = zeros(1,860); % NOTE: 860 is # of lagrangian pts being analyzed in tube
fTanAvg = fMagAvg;
fNormAvg= fMagAvg;

ct = 0;
for i=(start+dump_int):dump_int:finish
        
        iP = i - dump_int; % previous time-step index variable
    
        ct = ct + 1;       % Counter for storage indexing
        
        fprintf('\nAnalyzing %d of %d\n',i/dump_int,finish/dump_int)
    
        % Points to desired data viz_IB2d data file for CURRENT time-step
        if i<10
           numSim = ['000', num2str(i) ];
        elseif i<100
           numSim = ['00', num2str(i) ];
        elseif i<1000
           numSim = ['0', num2str(i)];
        else
           numSim = num2str(i);
        end
        
        % Points to desired data viz_IB2d data file for PREVIOUS time-step data
        if iP<10
           numSim_Prev = ['000', num2str(iP) ];
        elseif iP<100
           numSim_Prev = ['00', num2str(iP) ];
        elseif iP<1000
           numSim_Prev = ['0', num2str(iP)];
        else
           numSim_Prev = num2str(iP);
        end
        
        % Imports immersed boundary positions %
        [xLag,yLag] = give_Lag_Positions(pathForce,numSim);
        [xLag_Prev,yLag_Prev] = give_Lag_Positions(pathForce,numSim_Prev);

        % Computes spatial distance between current and previous Lag. Pts.
        dist_Vec = ( (xLag - xLag_Prev ).^2 + ( yLag - yLag_Prev ).^2 ).^(1/2);

        % Computes magnitude of velocity between current and previous Lag. Pts.
        speed_Vec = dist_Vec / dt;
        

        % Imports Lagrangian Pt. FORCE (magnitude) DATA %
        %                      DEFINITIONS 
        %
        %      fX_Lag: forces in x-direction on boundary
        %      fY_Lag: forces in y-direction on boundary
        %       fLagMag: magnitude of force at boundary
        %   fLagNorm: magnitude of NORMAL force at boundary
        %   fLagTan: magnitude of TANGENT force at boundary
        %
        [fLagMag,fLagNorm,fLagTan] = import_Lagrangian_Force_Data(pathForce,numSim);
        
        % Compute Cost of Transport at Particular Time-Step ( COT_i = |V_i||F_i| / L )
        COT_Vec(ct) = sqrt( speed_Vec'*speed_Vec )*sqrt( fLagMag'*fLagMag ) / L; 
        
        % Compute work for each Lag. Point
        Work_Vec = dist_Vec .* fLagMag;
        
        % Compute total work over time-step
        Work(ct) = sum(Work_Vec);
 
end

% Compute Cost of Transport Quantity for whole Simulation
COT_Sim = mean( COT_Vec );           % Computes average COT across simulation

% Compute Average Work Done in Simulation
Work_Sim = mean(Work);               % Computes average work

% PLOT WORK / COT over time!
if plot_Flag == 1

    % Plot the WORK done over time!
    figure(1)
    lw = 4;
    ms = 20;
    fs = 20;
    plot(start+1:dump_int:finish,Work,'-.','MarkerSize',ms,'LineWidth',lw); 
    ylabel('Work: SUM (fMag)_i*(dist)_i');
    xlabel('time-step');
    set(gca,'FontSize',fs);

    % Plot the COST OF TRANSPORT over time!
    figure(2)
    lw = 4;
    ms = 20;
    fs = 20;
    plot(start+1:dump_int:finish,COT_Vec,'-.','MarkerSize',ms,'LineWidth',lw); 
    ylabel('COT = |V_i||F_i|');
    xlabel('time-step');
    set(gca,'FontSize',fs);

end

%strName = [sim '_Work'];
%print_Matrix_To_Txt_File(fNormAvg,strName)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    ENDS MAIN PART OF CODE !!!                         %
%             BELOW IS HOW IBAMR DATA IS READ IN / ANALYZED             %
%            BELOW THAT ARE ADDITIONAL SCRIPTS FOR ANALYSIS             %
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: function to analyze FORCE (hier_data_IB2d) data from IBAMR
%           simulations / print them to be opened in VisIT
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_IBAMR_hier_data_to_VTK(path_name,starting_time,time_step,final_time,timestep,folder_Name)


% FOR THIS CODE YOU MUST TELL IT:
%   1. Path to desired data (line 29)
%   2. Starting index of desired data, e.g., F.01000 -> starting_time = 1000 
%   3. Ending index of desired data to analyze, e.g, F.28000 -> final_time = 28000 
%   4. The interval between successive data dumps, e.g., if F.01000 and next is F.02000, and next is F.03000 -> time_step = 1000 
%   5. Name of folder for new tangential, magnitude of force data 
%   5. It should be automated otherwise.

% WHAT IT OUTPUTS: 
%   1. Dumps new force data into the folder you named in Line 59. 
%   2. Run script save_Work_Cost_Of_Transport.m to compute COT for particular simulation


% Saving indices (to only grab subset of total pts!)
top_inds = [616:1:676 0:1:307 985:1:1045];
bot_inds = [1046:1:1106 308:1:615 1415:1:1475];

% To go from IBAMR indexing (starts at 0) to IB2d indexing starts at 1
top_inds = top_inds + 1;
bot_inds = bot_inds + 1;

% For Re-dimensionalization (characteristic values) *if necessary* %
density= 1025;   % fluid dynamic density (from IBAMR)
velocity=0.0075; % characteristic velocity (m/s^2)
length=  208e-6; % characteristic length (m)

% Change to Folder To Store .vtk Files
mkdir([path_name,folder_Name]);

% Print program information to screen
fprintf('\n\n              <--*** PRINTS hier_IB2d_data to VTK format ***--> \n\n');

% Loop Over All The Simulation Data
for k=starting_time:time_step:final_time

    fprintf('->Printing: %d of %d\n',k,final_time);
    
    % Points to FileID for Specific Timestep
    if k < 10000
        fileID = fopen([path_name 'F.0' num2str(k)])  %This line opens the file of interest and gives it an FileID #
        fileID2 = fopen([path_name 'X.0' num2str(k)])
    else
        fileID = fopen([path_name 'F.' num2str(k)]);  %This line opens the file of interest and gives it an FileID #
        fileID2 = fopen([path_name 'X.' num2str(k)]);
    end
    
    % Read in Force Data
    [F_Lag,~] = read_in_IBAMR_hier_data(fileID);
       
    % Read in Lagrangian Positions
    [lagPts,~] = read_in_IBAMR_hier_data(fileID2);

    % Give me desired pts for analysis! (selecting to compute forces over peristaltic region of points only)
    [lagPts,F_Lag,Npts] = give_Me_Desired_Subset_of_Pts(lagPts,F_Lag,top_inds,bot_inds);
    
    % Compute Normal/Tangential Vectors
    [nX,nY,sqrtN] = give_Me_Lagrangian_Normal_Vectors(Npts,lagPts);
    [tX,tY,~] = give_Me_Lagrangian_Tangent_Vectors(Npts,nX,nY,sqrtN);
    
    % Project Force Data onto Normal / Tangent Vectors
    [F_Tan,F_Normal] = give_Tangent_and_Normal_Force_Projections(Npts,F_Lag,nX,nY,tX,tY);
    
    % Re-dimensionalize the Forces (since the force is the non-dimensional force, or force coefficient)
    %F_Tan = -2*F_Tan/(density*velocity^2*length);
    %F_Normal = -2*F_Normal/(density*velocity^2*length);
    
    % Compute Colormap Force Magnitude Scalings
    cutoff_Yes = 0; % Don't cutoff top and bottom (use for saturation, if desired)
    [F_Tan_Mag,F_Normal_Mag] = give_Force_Magnitude_Scalings(Npts,F_Tan,F_Normal,cutoff_Yes);
    
    % Print Force Data to VTK Format
    print_Forces_to_VTK([path_name,folder_Name],num2str(k),lagPts,F_Lag,F_Tan_Mag,F_Normal_Mag)
    
    % Save (x,y) Positions to VTK Format
    print_XY_Positions_to_VTK([path_name,folder_Name],num2str(k),lagPts);
    
    % Plots stuff in MATLAB (previously used for testing only)
    print_Lagrangian_Mesh_and_Forces(lagPts,F_Tan_Mag,F_Normal_Mag);

end

fprintf('\n\nDone Converting to .vtk format. Happy computing and enjoy your day!\n\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    ENDS MAIN PART OF CODE !!!                         %
%             BELOW ARE SCRIPTS TO COMPUTE DESIRED FORCES               %
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Give me desired pts for analysis!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lagPtsNew,F_LagNew,Npts] = give_Me_Desired_Subset_of_Pts(lagPts,F_Lag,top_inds,bot_inds)

% Grabs Full Row in each Matrix for a Specific Index corresponding to the TOP of tube
ct = 0;
for j=1:length(top_inds);
   
   ct = ct + 1;       % Saving info counter
   
   ind = top_inds(j); % Chosen index along top of tube (according to IBAMR)
 
   % Grabbing corresponding information from IBAMR and reordering
   lagPtsNew(ct,:) = lagPts( ind , :);
   F_LagNew(ct,:) = F_Lag( ind, :);
   
end

% Grabs Full Row in each Matrix for a Specific Index corresponding to the BOTTOM of tube
for j=1:length(bot_inds);
   
   ct = ct + 1;       % Saving info counter
   
   ind = bot_inds(j); % Chosen index along top of tube (according to IBAMR)
 
   % Grabbing corresponding information from IBAMR and reordering
   lagPtsNew(ct,:) = lagPts( ind , :);
   F_LagNew(ct,:) = F_Lag( ind, :);
   
end

Npts = length(lagPtsNew(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: print (X,Y) coordinates of Lag. Data to .VTK Format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_XY_Positions_to_VTK(folder_name,strNUM,lagPts)

cd(folder_name); %change directory to force_VTK_Files folder

% Define name for 'filename'.vtk 
filename = ['lagsPts.' strNUM '.vtk'];

% Save (x,y) pts to .vtk format
savevtk_points( lagPts, filename, 'lagPts')

cd .. % Get out of force_VTK_Files


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints matrix vector data to vtk formated file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_points( X, filename, vectorName)

%X is matrix of size Nx3

N = length( X(:,1) );


%TRY PRINTING THEM AS UNSTRUCTURED_GRID
file = fopen (filename, 'w');
fprintf(file, '# vtk DataFile Version 2.0\n');
fprintf(file, [vectorName '\n']);
fprintf(file, 'ASCII\n');
fprintf(file, 'DATASET UNSTRUCTURED_GRID\n\n');
%
fprintf(file, 'POINTS %i float\n', N);
for i=1:N
    fprintf(file, '%.15e %.15e %.15e\n', X(i,1),X(i,2),0);
end
fprintf(file,'\n');
%
fprintf(file,'CELLS %i %i\n',N,2*N); %First: # of "Cells", Second: Total # of info inputed following
for s=0:N-1
    fprintf(file,'%i %i\n',1,s);
end
fprintf(file,'\n');
%
fprintf(file,'CELL_TYPES %i\n',N); % N = # of "Cells"
for i=1:N
   fprintf(file,'1 '); 
end
fprintf(file,'\n');
fclose(file);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: print force data to VTK format
%           -> each lag. pt (x,y) has associated scalar force with it
%           -> force magnitude, normal force, tangential force
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Forces_to_VTK(folder_Name,strNUM,lagPts,F_Lag,F_Tan_Mag,F_Normal_Mag)


cd(folder_Name); %change directory to force_VTK_Files folder

    fMagName = ['fMag.' strNUM '.vtk'];
    fNormalName = ['fNorm.' strNUM '.vtk'];
    fTangentName = ['fTan.' strNUM '.vtk'];

    fLagMag = sqrt( F_Lag(:,1).^2 + F_Lag(:,2).^2 ); % Compute magnitude of forces on boundary

    savevtk_points_with_scalar_data( lagPts, fLagMag, fMagName, 'fMag');
    savevtk_points_with_scalar_data( lagPts, F_Normal_Mag, fNormalName, 'fNorm');
    savevtk_points_with_scalar_data( lagPts, F_Tan_Mag, fTangentName, 'fTan');

cd .. % Get out of force_VTK_Files


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints Lagrangian pt data w/ associated scalar to vtk formated file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_points_with_scalar_data( X, scalarArray, filename, vectorName)

%X is matrix of size Nx3

N = length( X(:,1) );


%TRY PRINTING THEM AS UNSTRUCTURED_GRID
file = fopen (filename, 'w');
fprintf(file, '# vtk DataFile Version 2.0\n');
fprintf(file, [vectorName '\n']);
fprintf(file, 'ASCII\n');
fprintf(file, 'DATASET UNSTRUCTURED_GRID\n\n');
%
fprintf(file, 'POINTS %i float\n', N);
for i=1:N
    fprintf(file, '%.15e %.15e %.15e\n', X(i,1),X(i,2),0);
end
fprintf(file,'\n');
%
fprintf(file, 'POINT_DATA   %d\n', N);
fprintf(file, ['SCALARS ' vectorName ' double\n']);
fprintf(file, 'LOOKUP_TABLE default\n');
fprintf(file, '\n');
    for i=1:N
        fprintf(file, '%d ', scalarArray(i,1));
        fprintf(file, '\n');
    end

fclose(file);
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: read in data from hier_data from IBAMR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F_Lag,Npts] = read_in_IBAMR_hier_data(fileID)

%Scans the file for numbers (skipping 'Processor' lines) and puts them in a cell called C
C=textscan(fileID,'%f','HeaderLines',3,'delimiter',' ','commentstyle','P'); 

Forces=C{1}; %Converts cell C to matrix A

F_Lag(:,1) = Forces(1:2:end); % Converts data to x-directed Forces
F_Lag(:,2) = Forces(2:2:end); % Converts data to y-directed Forces
Npts = length(Forces(1:2:end));    % # of Lagrangian Pts.

fclose(fileID); %close the FileID #


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes Lagrangian Derivatives for Normal/Tangential Vector Computation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xL_s,yL_s] = give_Me_Lagrangian_Derivatives(ds,Npts,X,Y)

xL = X;
yL = Y;

xL_s = zeros(Npts,1);
yL_s = zeros(Npts,1);

% NOTE: dS = true distance between end pts in symmetric 3-pt stencil, so no
%            need for 2*ds in denominator.

for i=1:Npts
    if i==1
       dS = sqrt( ( xL(2) - xL(end) )^2 + ( yL(2) - yL(end) )^2 );
       xL_s(1) = ( xL(2) - xL(end) ) / (dS); 
       yL_s(1) = ( yL(2) - yL(end) ) / (dS);
    elseif i<Npts
       ds = sqrt( ( xL(i+1) - xL(i-1) )^2 + ( yL(i+1) - yL(i-1) )^2 );
       xL_s(i) = ( xL(i+1) - xL(i-1) ) / (dS); 
       yL_s(i) = ( yL(i+1) - yL(i-1) ) / (dS);
    else
       dS = sqrt( ( xL(1) - xL(end-1) )^2 + ( yL(1) - yL(end-1) )^2 );
       xL_s(i) = ( xL(1) - xL(end-1) ) / (dS); 
       yL_s(i) = ( yL(1) - yL(end-1) ) / (dS);
    end
    
    % TROUBLE! If both derivatives are zero! Default to previous value.
    if ( xL_s(i) == 0 ) && ( yL_s(i) == 0 )
       xL_s(i) = xL(i-1); 
       yL_s(i) = yL(i-1);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes Lagrangian UNIT Normal Vectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nX,nY,sqrtN] = give_Me_Lagrangian_Normal_Vectors(Npts,lagPts)

% Lagrangian Pts
X = lagPts(:,1); % x-Lagrangian Values
Y = lagPts(:,2); % y-Lagrangian Values

% Compute Lagrangian Spacing
%ds = sqrt( ( X(3)-X(4) )^2 + ( Y(3)-Y(4) )^2 );
ds = 1; % dummy

% Gives Lagrangian Derivatives
[xL_s,yL_s] = give_Me_Lagrangian_Derivatives(ds,Npts,X,Y);

sqrtN = sqrt( (xL_s).^2 + (yL_s).^2 );

nX = ( yL_s ) ./ sqrtN;
nY = ( -xL_s) ./ sqrtN;

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes Lagrangian UNIT Tangent Vectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tX,tY,sqrtN] = give_Me_Lagrangian_Tangent_Vectors(Npts,nX,nY,sqrtN)

% Allocate storage
tX = zeros( size(nX) );
tY = zeros( size(nY) );

% Rotate normal vectors to get tangent vectors
ang = -pi/2; % Rotate CW by 90 degrees
for i=1:Npts
    tX(i) = nX(i)*cos(ang) - nY(i)*sin(ang);
    tY(i) = nX(i)*sin(ang) + nY(i)*cos(ang);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes force vector projections onto the tangent and normal
%           vectors! 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [F_Tan,F_Normal] = give_Tangent_and_Normal_Force_Projections(Npts,F_Lag,nX,nY,tX,tY)

% Allocate Storage
F_Tan = zeros(Npts,2); F_Normal = F_Tan;

Fx = F_Lag(:,1); % Forces in x-direction
Fy = F_Lag(:,2); % Forces in y-direction

for i=1:Npts

    % Compute dot product between force vector and tangent vector
    tanVec_dotProd = ( Fx(i)*tX(i) + Fy(i)*tY(i) ) / sqrt( tX(i)*tX(i) + tY(i)*tY(i) );
    F_Tan(i,1) = tanVec_dotProd * ( tX(i) );
    F_Tan(i,2) = tanVec_dotProd * ( tY(i) );
    
    % Compute dot product between force vector and normal vector
    normalVec_dotProd = ( Fx(i)*nX(i) + Fy(i)*nY(i) ) / sqrt( nX(i)*nX(i) + nY(i)*nY(i) );
    F_Normal(i,1) = normalVec_dotProd * ( nX(i) );
    F_Normal(i,2) = normalVec_dotProd * ( nY(i) );
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: scales the force matrices by desired percentiles of each in magnitude 
%           for colormap scalings 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MagTan,MagNormal] = give_Force_Magnitude_Scalings(Npts,F_Tan,F_Normal,cutoff_Yes)

% Allocate Storage
MagTan = zeros( length( F_Tan(:,1) ),1 );
MagNormal = MagTan;

% Find the magnitude of the force in each 
for i=1:Npts
    MagTan(i,1) = sqrt( F_Tan(i,1)^2  + F_Tan(i,2)^2 );
    MagNormal(i,1)= sqrt( F_Normal(i,1)^2 + F_Normal(i,2)^2);
end

if cutoff_Yes == 1
    
    % Finds Percentiles for Forces for Cutoff Pts.
    prc90_T = prctile(MagTan,90);
    prc90_N = prctile(MagNormal,90);
    prc10_T = prctile(MagTan,10);
    prc10_N = prctile(MagNormal,10);

    % "Cutoff Threshold for the forces" via if-elseif statements by desired percentiles. 
    for i=1:Npts

        mT = MagTan(i);
        mN = MagNormal(i);

        if mT >= prc90_T
            MagTan(i) = prc10_T;
        elseif mT <= prc10_T
            MagTan(i) = prc10_T;
        end

        if mN >= prc90_N
            MagNormal(i) = prc10_N;
        elseif mN <= prc10_N
            MagNormal(i) = prc10_N;
        end

    end
    
end % Ends scaling

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: gives (x,y) positions of the immersed boundary at a single step
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = give_Lag_Positions(path,numSim)

analysis_path = pwd;

[xLag,yLag] = read_Lagrangian_Data_From_vtk(path,numSim);

cd(analysis_path);

clear analysis_path;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: imports all Eulerian Data at a single step
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function [fX_Lag,fY_Lag,fLagMag,fLagNorm,fLagTan] = import_Lagrangian_Force_Data(path,numSim)
function [fLagMag,fLagNorm,fLagTan] = import_Lagrangian_Force_Data(path,numSim)


% read in Mag. of Force %
%strChoice = 'fX_Lag'; 
%fX_Lag =  read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Tangential Force %
%strChoice = 'fY_Lag'; 
%fY_Lag = read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Force %
strChoice = 'fMag'; 
fLagMag =  read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Tangential Force %
strChoice = 'fTan'; 
fLagNorm = read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Normal Force %
strChoice = 'fNorm';
fLagTan =  read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);
 

clear strChoice first;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in (x,y) positions of the immersed boundary from .vtk
%           format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Fdata = read_Force_Scalar_Data_From_vtk(path,simNums,strChoice)

cd(path);

filename = [strChoice '.' num2str(simNums) '.vtk'];  % desired lagPts.xxxx.vtk file

fileID = fopen(filename);
if ( fileID== -1 )
    error('\nCANNOT OPEN THE FILE!');
end

str = fgets(fileID); %-1 if eof
if ~strcmp( str(3:5),'vtk');
    error('\nNot in proper VTK format');
end

% read in the header info %
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);

% stores # of Lagrangian Pts. as stated in .vtk file
numLagPts = sscanf(str,'%*s %f %*s',1); 

for i=1:numLagPts+5 % +5 to get to FORCE DATA
    str = fgets(fileID);
end

% read in the vertices %
[mat,count] = fscanf(fileID,'%f',numLagPts);
if count ~= numLagPts
   error('\nProblem reading in Lagrangian Pts. Force Data'); 
end

mat = reshape(mat, 1, count/1); % Reshape vector -> matrix (every 3 entries in vector make into matrix row)
forces = mat';              % Store vertices in new matrix

fclose(fileID);               % Closes the data file.

Fdata = forces(:,1);         % magnitude of the force

cd ..;                        % Change directory back to ../hier_IB2d_data/ directory

clear mat str filename fileID count;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints the Lagrangian Mesh with Scalar Value in MATLAB 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Mesh_and_Forces(lagPts,F_Tan_Mag,F_Normal_Mag)

    subplot(1,2,1)
    %interval = [250:7:1500 2600:7:3400]; % for testing desired interval from original .vertex IBAMR indexing
    interval = 1:length( lagPts(:,1) );   % gives new length after extracting only relevant indices (those in peristaltic region)
    scatter(lagPts(interval,1),lagPts(interval,2),55,F_Tan_Mag(interval),'filled');
    title('F_{Tangent}');
    subplot(1,2,2)
    scatter(lagPts(interval,1),lagPts(interval,2),55,F_Normal_Mag(interval),'filled');
    title('F_{Normal}')
    pause(0.5);
    clf;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION prints matrix to file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Matrix_To_Txt_File(a,strName)

nameTxt = [strName '.txt'];

fid = fopen(nameTxt, 'wt'); % Open for writing
for i=1:size(a,1)
   fprintf(fid, '%d ', a(i,:));
   fprintf(fid, '\n');
end
fclose(fid);

%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled  
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nick.battista@unc.edu
% Date Created: May 27th, 2015
% Institution: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in (x,y) positions of the immersed boundary from .vtk
%           format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xLag,yLag] = read_Lagrangian_Data_From_vtk(path,simNums)


cd(path);

filename = ['lagsPts.' num2str(simNums) '.vtk'];  % desired lagPts.xxxx.vtk file

fileID = fopen(filename);
if ( fileID== -1 )
    error('\nCANNOT OPEN THE FILE!');
end

str = fgets(fileID); %-1 if eof
if ~strcmp( str(3:5),'vtk');
    error('\nNot in proper VTK format');
end

% read in the header info %
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);

% stores # of Lagrangian Pts. as stated in .vtk file
numLagPts = sscanf(str,'%*s %f %*s',1); 

% read in the vertices %
[mat,count] = fscanf(fileID,'%f %f %f',3*numLagPts);
if count ~= 3*numLagPts
   error('\nProblem reading in Lagrangian Pts.'); 
end

mat = reshape(mat, 3, count/3); % Reshape vector -> matrix (every 3 entries in vector make into matrix row)
vertices = mat';              % Store vertices in new matrix

fclose(fileID);               % Closes the data file.

xLag = vertices(:,1);         % x-Lagrangian Pts.
yLag = vertices(:,2);         % y-Lagrangian Pts.

cd ..;                        % Change directory back to ../viz_IB2d/ directory

clear mat str filename fileID count;

