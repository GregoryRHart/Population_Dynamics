function Tcell_traj(directory)
% This function reads in the data from the T cells and plot the
% trajectories. 

tic;
%% parameters
fsize=12;
if nargin < 1
    directory = '';
end
if(~isempty(directory) && directory(end) ~= filesep)
    directory=[directory filesep];
end

% sizing screen for figure generation
set(0,'Units','pixels')
scnsize = [1,1,1920,1080];
pos_E = [2*scnsize(3)/3, scnsize(4)/3, scnsize(3)/3, scnsize(4)/3];

% resIdx for model
resIdx=cell(1,1);
m=0;                                    % at end of file read, m = # sites
mutidx=zeros(1,1);
fin_resIdx = fopen([directory 'resIdx.dat'],'rt');
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

traj_file=[directory 'Tcell_traj.dat'];

%% Reading in J data

% get file size array for preallocation
finfo = dir(traj_file);
fid = fopen(traj_file,'r');
intsize = fread(fid,1,'int32')*8;
floatsize = fread(fid,1,'int32')*8;
n_epitope = fread(fid,1,'int32');
n_snap = (finfo.bytes-12)/((intsize + floatsize*n_epitope*3+floatsize*n_epitope)/8);


% reading array
t = zeros(n_snap,1);
traj_array = zeros(n_snap,n_epitope*4);

for snap=1:n_snap
    t(snap) = fread(fid, 1, ['int' num2str(intsize)]);
    for i=1:n_epitope
        traj_array(snap,((i-1)*4+1):i*4-1) = fread(fid, 3, ['float' num2str(floatsize)]);
        traj_array(snap,i*4) = fread(fid, 1, ['float' num2str(floatsize)]);
    end
end

fclose(fid);

%% Making figures

% plotting trajectories
fig_handle = figure('numbertitle','off','name','T cell traj','Position',pos_E);

plot(t,traj_array(:,1:4:end)/10,':');
hold on;
plot(t,traj_array(:,2:4:end));
plot(t,traj_array(:,3:4:end),'--');
xlabel('iteration ','fontsize',fsize);
ylabel('# T cells','fontsize',fsize);
set(gca,'fontsize',fsize);
labels = {};
for i=1:n_epitope
     labels = [labels, {['1/10 naive_' num2str(i)]}];
end
for i=1:n_epitope
     labels = [labels, {['effector_' num2str(i)]}];
end
for i=1:n_epitope
     labels = [labels, {['memory_' num2str(i)]}];
end
legend(labels);

saveas(gcf,traj_file(1:end-4),'fig');
print(gcf,'-djpeg',[traj_file(1:end-4),'.jpg']);

close(fig_handle);

fig_handle = figure('numbertitle','off','name','immunogen traj','Position',pos_E);

plot(t,traj_array(:,4:4:end));
xlabel('iteration ','fontsize',fsize);
ylabel('chi*I','fontsize',fsize);
set(gca,'fontsize',fsize);
labels = {};
for i=1:n_epitope
     labels = [labels, {['sum chi*I_' num2str(i)]}];
end
legend(labels);

saveas(gcf,[traj_file(1:end-14) 'immunogen_traj'],'fig');
print(gcf,'-djpeg',[traj_file(1:end-14),'immunogen_traj.jpg']);

close(fig_handle);

toc;
end
