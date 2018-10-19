function [n_pop, n_snap, t, traj_array] = load_MCseqs(m, filepath, filename)
% This function reads in data from the file 'filename' found in the
% directory 'filepath'. 'filename' could be one of a number of files
% producted by the Population Dynamics code. The is a file for the whole
% protein (MC_seqs.dat) and files for each of the epitopes targeted in the
% run. They are binary files that contains the time step and the sequences
% of the population (full protein or a specific epitope) from that time
% step. The first thing in the file are 2 32-bit integers that give the
% size of an integer and a short integer, then there is an integer with
% that contains the population size.
%
% VAR IN
% m         - Length of the protein or epitope in the file.
% filepath  - filepath is an optional input arguement that allows the user
%             to specify the directory they would like to read
%             'filename' from. If it is not given the function reads 
%             from the current working directory.
% filename  - The name of the file to read in. This should be either
%             MC_seqs.dat or an epitope file which will be named with the
%             position, epitope sequence, and restricting HLA. If it is not
%             given MC_seqs.dat is assumed.
%
% VAR OUT
% n_pop        - This is the size of the population recorded in the file.
% n_snap       - This is the number of time steps recorded in the file. It
%                is calculated from the size of the file, the size of
%                integers and floats.
% t            - t is an array of integers storing the steps at which data
%                was recorded.
% traj_array   - traj_array is an array of short ints that n_snap by n_pop
%                by m in size. It stores the population sequences at each 
%                time step.
%
% FILE IN
% 'filename'   - This script reads in 'filename'. 'filename is one of the 
%                files output by the Population Dynamics code.
%
% FILE OUT
% none
%
% Written by:
% GRH 22 Feb 2017
%

%% Setup

tic;
% Default filename is MC_seqs.dat
% Default filepath is the current directory
if nargin < 3
    filename = 'MC_seqs.dat';
    if nargin < 2
        filepath ='';
    end
end
% Make sure that the filepath end in a file seperator
if(~isempty(filepath) && filepath(end) ~= filesep)
    filepath=[filepath filesep];
end

%% Reading in MC_seqs gilsp data

traj_file=[filepath filename];

% get file size array for preallocation
finfo = dir(traj_file);
fid = fopen(traj_file,'r');
cyclesize = fread(fid,1,'int32')*8;
seqsize = fread(fid,1,'int32')*8;
n_pop = fread(fid,1,['int' num2str(cyclesize)]);
n_snap = (finfo.bytes-8-cyclesize/8)/((cyclesize + seqsize*n_pop*m)/8);


% pre-allocate arrays
t = zeros(n_snap,1);
traj_array = zeros(n_snap,n_pop,m);

% reading array
for snap=1:n_snap
    t(snap) = fread(fid, 1, ['int' num2str(cyclesize)]);
    traj_array(snap,:,:) = reshape(fread(fid, m*n_pop, ['int' num2str(seqsize)]),[m,n_pop])';
end

fclose(fid);

disp(['Load MC_seqs run time: ', num2str(toc)])