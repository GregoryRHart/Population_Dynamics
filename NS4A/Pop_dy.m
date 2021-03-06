function Pop_dy(filepath,P2traj,T)
% This function reads in the data from the popualtion dynamics code and
% calculates P1, P2, and entropy and makes the corresponding plots.

tic;
%% parameters
if nargin < 3
    T = 1;
if nargin < 2
    P2traj = false;
    if nargin < 1
        filepath ='';
    end
end
end
if(~isempty(filepath) && filepath(end) ~= filesep)
    filepath=[filepath filesep];
end

fsize=12;

% sizing screen for figure generation
set(0,'Units','pixels')
scnsize = [1,1,1920,1080];
pos_N = [scnsize(3)/3, 2*scnsize(4)/3, scnsize(3)/3, scnsize(4)/3];
pos_S = [scnsize(3)/3, 5, scnsize(3)/3, scnsize(4)/3];


%% read in resIdx

% resIdx from old model
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

n_P1 = sum(nRes);
n_P2 = 0;
for i=1:m-1
   n_P2 = nRes(i)*sum(nRes(i+1:end));
end
%% Reading in MC_seqs data

traj_file=[filepath 'MC_seqs.dat'];

% get file size array for preallocation
finfo = dir(traj_file);
fid = fopen(traj_file,'r');
cyclesize = fread(fid,1,'int32')*8;
seqsize = fread(fid,1,'int32')*8;
n_pop = fread(fid,1,['int' num2str(cyclesize)]);
n_snap = (finfo.bytes-8-cyclesize/8)/((cyclesize + seqsize*n_pop*m)/8);


% reading array
t = zeros(n_snap,1);
traj_array = zeros(n_snap,n_pop,m);


for snap=1:n_snap
    t(snap) = fread(fid, 1, ['int' num2str(cyclesize)]);
    traj_array(snap,:,:) = reshape(fread(fid, m*n_pop, ['int' num2str(seqsize)]),[m,n_pop])';
end

fclose(fid);

%% Calc P1 and P2 trajectories
P1_traj = zeros(n_snap,n_P1);
if(P2traj)
P2_traj = zeros(n_snap,n_P2);
else
    P2_traj = zeros(1,n_P2);
end

uniqueResIdx = 1:21;              % unique residue indices residing in MSA

for snap=1:n_snap
    % one-site residue ordering and target probabilities
    hist1=histc(squeeze(traj_array(snap,:,:)),uniqueResIdx);              % one-site mutational histograms
    [hist1,~]=sort(hist1,1,'descend');     % reordering the residue indices at each site in descending order of prevalence
    
    P1_model=cell(1,m);
    for i=1:m
        P1_model{i}=hist1(1:nRes(i),i)/n_pop;
    end
    
    
    % two-site target probabilities
    % \-> matrix indexing based on the one-site residue ordering in each pair
    if(P2traj || snap == n_snap)
    P2_model=cell(m,m);
    for i=1:m
        for j=i+1:m
            [edge_i,idx_sort_i]=sort(resIdx{i},'ascend');            % hist3 requires that edges must be monotonically increasing
            [edge_j,idx_sort_j]=sort(resIdx{j},'ascend');
            
            idx_unsort_i=nan(length(idx_sort_i),1);
            idx_unsort_j=nan(length(idx_sort_j),1);
            idx_unsort_i(idx_sort_i)=1:1:length(idx_sort_i);
            idx_unsort_j(idx_sort_j)=1:1:length(idx_sort_j);
            
            % hist3 freaks out if either dimension contains only one bin, this
            %  hack adds a dummy bin to the singleton dimension and then strips
            %  out the dummy col and/or row of the resulting histogram
            if ( length(edge_i)>1 )
                if ( length(edge_j)>1 )
                    p=hist3(squeeze(traj_array(snap,:,[i,j])),'Edges',{edge_i,edge_j});
                else % length(edge_j)==1
                    p=hist3(squeeze(traj_array(snap,:,[i,j])),'Edges',{edge_i,cat(1,edge_j,200184*edge_j)});
                    p=p(:,1);
                end
            else	% length(edge_i)==1
                if ( length(edge_j)>1 )
                    p=hist3(squeeze(traj_array(snap,:,[i,j])),'Edges',{cat(1,edge_i,170786*edge_i),edge_j});
                    p=p(1,:);
                else % length(edge_j)==1
                    p=hist3(squeeze(traj_array(snap,:,[i,j])),'Edges',{cat(1,edge_i,170786*edge_i),cat(1,edge_j,200184*edge_j)});
                    p=p(1,1);
                end
            end
            
            p=p(idx_unsort_i,idx_unsort_j)/n_pop;                     % reordering and normalizing
            P2_model{i,j}=p;
        end
    end
        if(P2traj)
            temp_idx = snap;
        else
            temp_idx = 1;
    end
    
    pos = 1;
    for i=1:m
        for j=i+1:m
            for q=1:nRes(j)
                    P2_traj(temp_idx,pos:(pos+(nRes(i))-1)) = P2_model{i,j}(:,q);
                pos = pos + nRes(i);
            end
        end
    end
        
    end
    
    pos = 1;
    for i=1:m
        P1_traj(snap,pos:(pos+nRes(i)-1)) = P1_model{i};
        pos = pos + nRes(i);
    end
    
end


%% Making figures

% plotting trajectories
fig_handle = figure('Position',pos_N);

plot(t,P1_traj);
xlabel('iteration ','fontsize',fsize);
ylabel('value ','fontsize',fsize);
set(gca,'fontsize',fsize);

saveas(gcf,[filepath 'P1_traj'],'fig');
print(gcf,'-djpeg',[filepath,'P1_traj.jpg']);

close(fig_handle);

if(P2traj)
% plotting trajectories
fig_handle = figure('Position',pos_S);

plot(t,P2_traj);
xlabel('iteration ','fontsize',fsize);
ylabel('value ','fontsize',fsize);
set(gca,'fontsize',fsize);

saveas(gcf,[filepath 'P2_traj'],'fig');
print(gcf,'-djpeg',[filepath,'P2_traj.jpg']);

close(fig_handle);
end

% P1_fit
traj_file=[filepath 'P1_target.dat'];
fid = fopen(traj_file,'rt');
P1_target = fgetl(fid);                 % reading line
P1_target = sscanf(P1_target,'%f');     % splitting line
fclose(fid);

fig_handle = figure('Position',pos_N);

scatter(P1_target,P1_traj(end,:),25,'r');
hold on
plot([0 1],[0 1],'-k');
hold off
xlabel('P1\_target ','fontsize',fsize);
ylabel('P1\_model (MC) ','fontsize',fsize);
set(gca,'fontsize',fsize); 
xlim([0 1.1*max([P1_target',P1_traj(end,:)])])
ylim([0 1.1*max([P1_target',P1_traj(end,:)])])

saveas(gcf,[filepath 'P1_fit'],'fig');
print(gcf,'-djpeg',[filepath 'P1_fit.jpg']);

close(fig_handle);





% P2_fit
traj_file=[filepath 'P2_target.dat'];
fid = fopen(traj_file,'rt');
P2_target = fgetl(fid);                 % reading line
P2_target = sscanf(P2_target,'%f');     % splitting line
fclose(fid);

fig_handle = figure('Position',pos_N);

scatter(P2_target,P2_traj(end,:),25,'r');
hold on
plot([0 1],[0 1],'-k');
hold off
xlabel('P2\_target ','fontsize',fsize);
ylabel('P2\_model (MC) ','fontsize',fsize);
set(gca,'fontsize',fsize); 
xlim([0 1.1*max([P2_target',P2_traj(end,:)])])
ylim([0 1.1*max([P2_target',P2_traj(end,:)])])

saveas(gcf,[filepath 'P2_fit'],'fig');
print(gcf,'-djpeg',[filepath 'P2_fit.jpg']);

close(fig_handle);

clear P1_traj
clear P2_traj

%% Calculate Entropy

entropy_traj = zeros(length(t),3);%10);
parfor snap=1:n_snap
    P = 0;
    seq = squeeze(traj_array(snap,:,:));
    idx = randperm(n_pop);
    for i=1:3%10
        sample = sortrows(seq(idx(ceil(1:(n_pop*(i*.1 + .7)))),:));
        temp = sample(1,:);
        for j=ceil(1:(n_pop*(i*.1 + .7)))
            if(all(sample(j,:)==temp))
               P = P + 1;
            else
                temp = sample(j,:);
                P = P/length(sample);
            if(P~=0)
                entropy_traj(snap,i) = entropy_traj(snap,i) - P.*log(P);
            end
                P = 1;
            end
        end
    end
end

clear traj_array

%% Reading in fitness data

traj_file=[filepath 'pop_stats.dat'];

% reading array
fitness_array = zeros(n_snap,4);

fid = fopen(traj_file,'rt');
intsize = fread(fid,1,'int32')*8;
floatsize = fread(fid,1,'int32')*8;

for snap=1:n_snap    
    fread(fid, 1, ['int' num2str(intsize)]);
    fitness_array(snap,1) = fread(fid, 1, ['int' num2str(intsize)]);
    fitness_array(snap,2:end) = fread(fid, 3, ['float' num2str(floatsize)]);  
end

fclose(fid);
E = -T*log(fitness_array(:,2));
E_eff = -T*log(fitness_array(:,4));

plot(t,E - T*mean(entropy_traj,2),t,E_eff - T*mean(entropy_traj,2))

saveas(gcf,[filepath 'free_fit'],'fig');
print(gcf,'-djpeg',[filepath 'free_fit.jpg']);

S = mean(entropy_traj,2);
save([filepath 'freeFitness.mat'],'E','E_eff','S');


toc;
end
