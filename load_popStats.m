function [n_snap, t, fitness_array] = load_popStats(filepath)
% This function reads in data from the file pop_stats.dat found in the
% directory 'filepath'. pop_stats.dat is producted by the Population
% Dynamics code. It is a binary file that contains the time step, the
% population size, the average effective fitness of the progeny, the
% average fitness and effective fitness of the population, and the average
% energy and effective energy of the population. The first thing in the
% file are 2 32-bit integers that give the size of an integer and a float.
%
% VAR IN
% filepath  - filepath is an optional input arguement that allows the user
%             to specify the directory they would like to read
%             pop_stats.dat from. If it is not given the function reads 
%             from the current working directory.
%
% VAR OUT
% n_snap       - This is the number of time steps recorded in the file. It
%                is calculated from the size of the file, the size of
%                integers and floats.
% t            - t is an array of integers storing the steps at which data
%                was recorded.
% finess_array - fitness_array is an array of floats that n_snap by 6 in
%                size. It stores the population size, fitness and effective
%                fitness, and energy and effective energy at each time
%                step.
%
% FILE IN
% pop_stats.dat - This script reads in pop_stats.dat. pop_stats.dat is
%                  output by the Population Dynamics code.
%
% FILE OUT
% none
%
% Written by:
% GRH 21 Feb 2017
%

%% Setup

tic;
% Default filepath is the current directory
if nargin < 1
    filepath ='';
end
% Make sure that the filepath end in a file seperator
if(~isempty(filepath) && filepath(end) ~= filesep)
    filepath=[filepath filesep];
end

%% Load pop_stats.dat
traj_file=[filepath 'pop_stats.dat'];

% read in size information
finfo = dir(traj_file);
fid = fopen(traj_file,'rt');
intsize = fread(fid,1,'int32')*8;
floatsize = fread(fid,1,'int32')*8;
n_snap = (finfo.bytes-8)/((intsize*2+floatsize*5)/8);

% pre-allocate arrays
t = zeros(n_snap,1);
fitness_array = zeros(n_snap,6);

% read in fitness data
for snap=1:n_snap
    t(snap) = fread(fid, 1, ['int' num2str(intsize)]);
    fitness_array(snap,1) = fread(fid, 1, ['int' num2str(intsize)]);
    fitness_array(snap,2:end) = fread(fid, 5, ['float' num2str(floatsize)]);
end

fclose(fid);
disp(['Load pop_stats run time: ', num2str(toc)])