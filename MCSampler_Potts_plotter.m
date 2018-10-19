function MCSampler_Potts_plotter(filepath,P2traj)

tic;
% analysis of trajectories outputted by HReconstr_Potts_MATLAB to assess convergence 
if nargin < 2
    P2traj = false;
    if nargin < 1
        filepath ='';
    end
end
if(~isempty(filepath) && filepath(end) ~= filesep)
    filepath=[filepath filesep];
end
    
%% parameters
fsize=18;


% sizing screen for figure generation
set(0,'Units','pixels') 
scnsize = [1,1,1920,1080];
pos_N = [scnsize(3)/3, 2*scnsize(4)/3, scnsize(3)/3, scnsize(4)/3];
pos_S = [scnsize(3)/3, 5, scnsize(3)/3, scnsize(4)/3];

% resIdx for model
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

n_h = sum(nRes);
n_J = 0;
for i=1:m-1
   n_J = n_J + nRes(i)*sum(nRes(i+1:end));
end

traj_file=[filepath 'P1_model_traj.dat'];

% get file size array for preallocation
finfo = dir(traj_file);
fid = fopen(traj_file,'r');
intsize = fread(fid,1,'int32')*8;
floatsize = fread(fid,1,'int32')*8;
n_snap = (finfo.bytes-8)/((intsize + floatsize*n_h)/8);


% reading array
t = zeros(n_snap,1);
traj_array = zeros(n_snap,n_h);

for snap=1:n_snap
    t(snap) = fread(fid, 1, ['int' num2str(intsize)]);
    traj_array(snap,:) = fread(fid, n_h, ['float' num2str(floatsize)]);
end

fclose(fid);

% plotting trajectories
fig_handle = figure('Position',pos_N);

plot(t,traj_array);
xlabel('iteration ','fontsize',fsize);
ylabel('value ','fontsize',fsize);
set(gca,'fontsize',fsize);

saveas(gcf,[traj_file(1:end-4),'.fig']);
print(gcf,'-djpeg',[traj_file(1:end-4),'.jpg']);

close(fig_handle);



if(P2traj)
traj_file=[filepath 'P2_model_traj.dat'];

% get file size array for preallocation
finfo = dir(traj_file);
fid = fopen(traj_file,'r');
intsize = fread(fid,1,'int32')*8;
floatsize = fread(fid,1,'int32')*8;
n_snap = (finfo.bytes-8)/((intsize + floatsize*n_J)/8);


% reading array
t = zeros(n_snap,1);
traj_array = zeros(n_snap,n_J);

for snap=1:n_snap
    t(snap) = fread(fid, 1, ['int' num2str(intsize)]);
    traj_array(snap,:) = fread(fid, n_J, ['float' num2str(floatsize)]);
end

fclose(fid);

% plotting trajectories
fig_handle = figure('Position',pos_S);

plot(t,traj_array);
xlabel('iteration ','fontsize',fsize);
ylabel('value ','fontsize',fsize);
set(gca,'fontsize',fsize);

saveas(gcf,[traj_file(1:end-4),'.fig']);
print(gcf,'-djpeg',[traj_file(1:end-4),'.jpg']);

close(fig_handle);
end



% P1_fit
traj_file=[filepath 'P1_target.dat'];
fid = fopen(traj_file,'rt');
    P1_target = fgetl(fid);                 % reading line
    P1_target = sscanf(P1_target,'%f');     % splitting line
fclose(fid);

traj_file=[filepath 'P1_model.dat'];
fid = fopen(traj_file,'r');
floatsize = fread(fid,1,'int32')*8;
P1_model = fread(fid, n_h, ['float' num2str(floatsize)]);
fclose(fid);

fig_handle = figure('Position',pos_N);

scatter(P1_target,P1_model,25,'r');
hold on
plot([0 1],[0 1],'-k');
hold off
xlabel('P1\_target ','fontsize',fsize);
ylabel('P1\_model (MC) ','fontsize',fsize);
set(gca,'fontsize',fsize); 
xlim([0 1.1*max([P1_target',P1_model'])])
ylim([0 1.1*max([P1_target',P1_model'])])
text(.25*1.1*max([P1_target',P1_model']),.75*1.1*max([P1_target',P1_model']),...
             ['SSE = ' num2str(sum((P1_target-P1_model).^2))],'FontSize',fsize);

saveas(gcf,[filepath 'P1_fit_MC.fig']);
print(gcf,'-djpeg',[filepath 'P1_fit_MC.jpg']);
set(gcf,'PaperPositionMode','auto');
print(gcf,'-depsc2', [filepath 'P1_fit_MC.eps']);

close(fig_handle);





% P2_fit
traj_file=[filepath 'P2_target.dat'];
fid = fopen(traj_file,'rt');
    P2_target = fgetl(fid);                 % reading line
    P2_target = sscanf(P2_target,'%f');     % splitting line
fclose(fid);

traj_file=[filepath 'P2_model.dat'];
fid = fopen(traj_file,'r');
floatsize = fread(fid,1,'int32')*8;
P2_model = fread(fid, n_J, ['float' num2str(floatsize)]);
fclose(fid);

fig_handle = figure('Position',pos_N);

scatter(P2_target,P2_model,25,'r');
hold on
plot([0 1],[0 1],'-k');
hold off
xlabel('P2\_target ','fontsize',fsize);
ylabel('P2\_model (MC) ','fontsize',fsize);
set(gca,'fontsize',fsize); 
xlim([0 1.1*max([P2_target',P2_model'])])
ylim([0 1.1*max([P2_target',P2_model'])])

saveas(gcf,[filepath 'P2_fit_MC.fig']);
print(gcf,'-djpeg',[filepath 'P2_fit_MC.jpg']);
set(gcf,'PaperPositionMode','auto');
print(gcf,'-djpeg',[filepath 'P2_fit_MC.eps']);

close(fig_handle);
toc;
