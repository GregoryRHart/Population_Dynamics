function [m, resIdx, nRes, mutidx] = load_resIdx(filepath)
% This function reads in the file resIdx.dat found in the directory
% 'filepath'. The file resIdx.dat contains informaton for indexing the
% reduced proteins for the Potts fitting code. It is also used to translate
% between the actual amino acids and the fake amino acids used in may of
% our codes.
%
% VAR IN
% filepath  - filepath is an optional input arguement that allows the user
%             to specify the directory they would like to read resIdx.dat
%             from. If it is not given it would read from the current
%             working directory.
%
% VAR OUT
% m         - m is the number of mutating sites in the protein. It is also
%             the number of lines in resIdx.dat.
% resIdx    - resIdx is a cell array of size m. Each cell in resIdx
%             contains the amino acids that appear in that (fake) position.
%             The are listed has numbers according to MATLAB's Mapping
%             Amino Acid Letter Codes to Integers table.
% nRes      - nRes is an array of size m. Each element gives the number of
%             different amino acids that can appear in that (fake)
%             position.
% mutidx    - mutidx is an array of size m. Because resIdx.dat can contain
%             only protein positions that mutate, we end up with a fake
%             position index where conserved positions have been removed.
%             mutidx gives the real positions of these mutating sites. For
%             example mutidx(1) would give the position in the full protein
%             of the first non-conserved site.
%
% FILE IN
% resIdx.dat    - This script reads in resIdx.dat. resIdx.dat is created by
%                 P_calc_reg or P_calc_reg_weight.
%
% FILE OUT
% none
%
% Written by:
% GRH 16 Feb 2017
% Based on parts of existing scripts written by ALF June 2012
%

%% Load resIdx.dat

tic;
% Default filepath is the current directory
if nargin < 1
    filepath ='';
end
% Make sure that the filepath end in a file seperator
if(~isempty(filepath) && filepath(end) ~= filesep)
    filepath=[filepath filesep];
end

% resIdx & nRes
resIdx=cell(1,1);
m=0;                                    % at end of file read, m = # sites
mutidx=zeros(1,1);
fin_resIdx = fopen([filepath 'resIdx.dat'],'rt');
while ~feof(fin_resIdx)
    m=m+1;
    line = fgetl(fin_resIdx);       % reading line
    line = sscanf(line,'%f');   	% splitting line
    mutidx(m)=line(1);
    resIdx{1,m} = line(2:end)';
end
fclose(fin_resIdx);

nRes=nan(1,m);
for i=1:m
    nRes(i)=length(resIdx{i});
end

disp(['Load resIdx run time: ', num2str(toc)])