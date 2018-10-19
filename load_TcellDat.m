function [n_epitope, n_snap, t, traj_array] = load_TcellDat(filepath)
% This function reads in data from the file Tcell_traj.dat found in the
% directory 'filepath'. Tcell_traj.dat is producted by the Population
% Dynamics code. It is a binary file that contains the time step, number of
% naive, effector, and memory cells, and the average susceptibility. The
% first thing in the file are 3 32-bit integers that give the size of an
% integer, the size of a float, and the number of epitopes.
%
% VAR IN
% filepath  - filepath is an optional input arguement that allows the user
%             to specify the directory they would like to read 
%             Tcell_traj.dat from. If it is not given the function reads 
%             from the current working directory.
%
% VAR OUT
% n_epitope  - n_epitope is the number of epitopes that the file contains.
%              This is needed to pre-allocate the arrays.
% n_snap     - This is the number of time steps recorded in the file. It is
%              calculated from the size of the file, the number of epitopes
%              and the size of integers and floats.
% t          - t is an array of integers storing the steps at which data was
%              recorded.
% traj_array - traj_array is an array of floats that n_snap by 4*n_epitope
%              in size. It stores the number of naive cells, effector
%              cells, memory cells, and average susceptibility for each
%              epitope at each time step.
%
% FILE IN
% Tcell_traj.dat - This script reads in Tcell_traj.dat. Tcell_traj.dat is
%                  output by the Population Dynamics code. It will be empty
%                  if the code was run with no immune pressure.
%
% FILE OUT
% none
%
% Written by:
% GRH 20 Feb 2017
%

%% Load T-cell data output from the Populaton Dynamics code

tic;
% Default filepath is the current directory
if nargin < 1
    filepath ='';
end
% Make sure that the filepath end in a file seperator
if(~isempty(filepath) && filepath(end) ~= filesep)
    filepath=[filepath filesep];
end


traj_file=[filepath 'Tcell_traj.dat'];

% read in size information
finfo = dir(traj_file);
fid = fopen(traj_file,'r');
intsize = fread(fid,1,'int32')*8;
floatsize = fread(fid,1,'int32')*8;
n_epitope = fread(fid,1,'int32');
n_snap = (finfo.bytes-12)/((intsize + intsize*n_epitope*3+floatsize*n_epitope)/8);

% reading array
t = zeros(n_snap,1);
traj_array = zeros(n_snap,n_epitope*4);

for snap=1:n_snap
    % first element of every "line" is the time step
    t(snap) = fread(fid, 1, ['int' num2str(intsize)]);
    for i=1:n_epitope
        % Read in the number of each type of T-cell
        traj_array(snap,((i-1)*4+1):i*4-1) = fread(fid, 3, ['float' num2str(intsize)]);
        % Read in the average T-cell susceptibility
        traj_array(snap,i*4) = fread(fid, 1, ['float' num2str(floatsize)]);
    end
end

fclose(fid);
disp(['Load TcellDat run time: ', num2str(toc)])
