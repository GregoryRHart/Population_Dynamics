function plotGilProb(filepath,P2traj,T)
% This function data from many different runs of the population dynamics
% code and calculates P1, P2, and entropy and makes plots of the mean,
% and standard deviation.

%% Setup

tic;
if nargin < 3
    T = 1;
    if nargin < 2
        P2traj = false;
        if nargin < 1
            filepath ='.';
        end
    end
end
if(~isempty(filepath) && filepath(end) ~= filesep)
    filepath=[filepath filesep];
end

fsize=12;

maincolor = [0,0,1; 0,1,0; 1,0,0; 0,1,1; 1,0,1; 1,1,0; 0,0,0];
seccolor = maincolor.*.5+.5;

d = dir(filepath);
d = {d([d.isdir] == 1).name};
d = regexpi(d,'[0-9]*','match');
d = [d{:}];

%% Read in data

% read in ODE solution and get size information
[~, n_snap, t, fitness_array] = evalc('load_popStats([filepath ''ODE''])');
EODE = fitness_array(:,5);
E_effODE = fitness_array(:,6);
fODE = fitness_array(:,3);
f_effODE = fitness_array(:,4);

% pre-allocate data arrays
Eall = zeros(n_snap,length(d));
E_effall = zeros(n_snap,length(d));
fall = zeros(n_snap,length(d));
f_effall = zeros(n_snap,length(d));
%Sall = zeros(n_snap,length(d));

% read in Gillespie data
for i=1:length(d)
    [~, n_snap, ~, fitness_array] = evalc('load_popStats([filepath d{i}])');
    Eall(:,i) = fitness_array(:,5);
    E_effall(:,i) = fitness_array(:,6);
    fall(:,i) = fitness_array(:,3);
    f_effall(:,i) = fitness_array(:,4);
end


E = zeros(n_snap,1);
E_eff = zeros(n_snap,1);
f = zeros(n_snap,1);
f_eff = zeros(n_snap,1);
%S = zeros(n_snap,1);
for i=1:n_snap
    E(i) = mean(Eall(i,:));
    E_eff(i) = mean(E_effall(i,:));
    f(i) = mean(fall(i,:));
    f_eff(i) = mean(f_effall(i,:));
    %    S(i) = mean(Sall(i,:));
end

%% Plot fitness

t_steps = ceil((length(t)-9)/5);

t_scale = [];
f_scale = [];
load('/data/grhart2/Population_Dynamics/long/M003/NS5B/invitro.mat','t_scale','f_scale')
t = t*t_scale;
plot_measureTraj(Eall, E, EODE, 'Energy',0)
plot_measureTraj(E_effall, E_eff, E_effODE, 'Effective Energy', 0)
plot_measureTraj(fall*exp(f_scale), f*exp(f_scale), fODE*exp(f_scale), 'Fitness', 1)
plot_measureTraj(f_effall*exp(f_scale), f_eff*exp(f_scale), f_effODE*exp(f_scale), 'Effective Fitness', 1)
%plot_measureTraj(Sall, S, SODE, 'Entropy', 1)


%% read in resIdx

% resIdx for the model
[~, ~, resIdx, nRes, mutidx] = evalc('load_resIdx(filepath)');

%% Find sites of interest

% load epitopes from file
epitopes = importdata([filepath d{1} '/epitopes.dat']);
if(isstruct(epitopes))
    epitopes = epitopes.textdata;
end
start = zeros(length(epitopes),1);
finish = zeros(length(epitopes),1);
for i=1:length(epitopes)
    [~, epitopes{i}] = strtok(epitopes{i}, '_');
    [token, ~] = strtok(epitopes{i}, '_');
    [s, f] = strtok(token,'-');
    start(i) = str2num(s); %#ok<ST2NM>
    finish(i) = str2num(f(2:end)); %#ok<ST2NM>
end

% calculet the number of mutating sites in the epitope and the number of
% possible single and double mutatants.
n_site = zeros(length(epitopes),1);
n_P1 = zeros(length(epitopes),1);
n_P2 = zeros(length(epitopes),1);
for i=1:length(epitopes)
    n_site(i) = (finish(i)-start(i)) + 1;
    n_P1(i) = sum(nRes(start(i):finish(i)));
    for k=start(i):(finish(i)-1)
        n_P2(i) = n_P2(i) + nRes(k)*sum(nRes(k+1:end));
    end
end

%% Reading in MC_seqs gilsp data

for k=1:length(epitopes)
    disp(['epitope ' num2str(k)])
    % pre-allocate arrays to hold mutational requencies.
    P1_all = zeros(n_snap,n_P1(k),length(d));
    %P2_all = zeros(n_snap,n_P2(k),length(d));
    
    % loop through all runs
    for c_run=1:length(d)
        [~, n_pop, n_snap, ~, traj_array] = evalc('load_MCseqs((finish(k)-start(k)+1), [filepath d{c_run}], epitopes{k}(2:end))');
        [P1_traj, ~] = P_calc(n_snap, n_P1(k), n_P2(k), traj_array, resIdx, nRes, start(k), finish(k), n_pop, P2traj);
        P1_all(:,:,str2num(d{c_run})) = P1_traj(:,:,1); %#ok<ST2NM>
        %P2_all(:,:,str2num(d{c_run})) = P2_traj(:,:,1); %#ok<ST2NM>
    end
    
    % find mean behavoir
    P1 = mean(P1_all,3);
    %P2 = mean(P2_all,3);
    
    %% Reading in MC_seqs ODE data  
    
    [~, n_pop, n_snap, ~, traj_array] = evalc('load_MCseqs((finish(k)-start(k)+1), [filepath d{c_run}], epitopes{k}(2:end))');
    [P1_ODE, ~] = P_calc(n_snap, n_P1(k), n_P2(k), traj_array, resIdx, nRes, start(k), finish(k), n_pop, P2traj);
    
    %% Plot residue trajectories
    
    [~, epitopes{k}] = strtok(epitopes{k}, '_');
    [epi,~] = strtok(epitopes{k},'_');
    labels = [];
    ticks = [];
    pos = 1;
    
    f = figure('numbertitle','off','name',['P1 for epitope ' epi]);
    for i=1:n_site(k)
        for r=1:nRes(start(k)+(i-1))
            count = 1;
            Xdata = zeros(t_steps,100);
            Ydata = zeros(t_steps,100);
            for j=1:5:(length(t)-9)
                temp = reshape(P1_all(j:j+9,pos,:),[],1);
                [Ydata(count,:),Xdata(count,:)] =  hist(temp,100);
                Xdata(count,:) = Xdata(count,:)+(i-1)*2;
                Ydata(count,:) = (Ydata(count,:)/max(Ydata(count,:)))*10 + t(j);
                count = count + 1;
            end
            plot(Ydata',Xdata','color',maincolor(mod(r,size(maincolor,1))+1,:))
            hold on
            pos = pos + 1;
        end
        ticks = [ticks; (i-1)*2; (i-1)*2+.5; (i-1)*2 + 1]; %#ok<AGROW>
        labels = [labels; {'0'}; {[int2aa(resIdx{start(k)+(i-1)}(1)) num2str(mutidx(start(k)+(i-1)))]};{'1'}]; %#ok<AGROW>
    end
    pos = 1;
    for i=1:n_site(k)
        for r=1:nRes(start(k)+(i-1))
            plot(t,P1(:,pos)+(i-1)*2,'Color',maincolor(mod(r,size(maincolor,1))+1,:)*.7,'LineWidth',2)
            plot(t,P1_ODE(:,pos)+(i-1)*2,'--','Color',maincolor(mod(r,size(maincolor,1))+1,:),'LineWidth',2)
            pos = pos + 1;
        end
    end
    
    ylim([0, n_site(k)*2-1])
    set(gca, 'ytick', ticks, 'yticklabel', labels);
    
    xlabel('Weeks', 'fontsize',fsize)
    ylabel('frequency', 'fontsize',fsize)
    set(gca,'fontsize',fsize);
    saveas(f,[filepath 'P1_' epi],'fig');
    print(f,'-djpeg',[filepath 'P1_' epi '.jpg']);
    close(f);
    
    labels = [];
    ticks = [];
    f = figure('numbertitle','off','name',['P1 for epitope ' epi]);
    image(t,1:n_P1(k),P1','CDataMapping','scaled');
    colormap('gray');
    hold on;
    plot(t,zeros(length(t),1)+0.5,'LineWidth',4,'color',[.5,0,.5]);
    count = 1;
    for i=1:n_site(k)
        plot(t,zeros(length(t),1)+0.5 + sum(nRes(start(k):(start(k)+(i-1)))),'LineWidth',6,'color',[.5,0,.5]);
        for r=1:nRes(start(k)+(i-1))
            ticks = [ticks; count]; %#ok<AGROW>
            labels = [labels; {[int2aa(resIdx{start(k)+(i-1)}(r)) num2str(mutidx(start(k)+(i-1)))]}]; %#ok<AGROW>
            count = count + 1;
        end
    end
    xlabel('Weeks', 'fontsize',fsize)
    ylabel('frequency', 'fontsize',fsize)
    set(gca,'fontsize',fsize);
    set(gca, 'ytick', ticks, 'yticklabel', labels);
    saveas(f,[filepath 'P1_' epi '_bar'],'fig');
    print(f,'-djpeg',[filepath 'P1_' epi '_bar.jpg']);
    close(f);
end

toc;

    function plot_measureTraj(all_traj, mean_traj, ODE_traj, plot_type, log_flag)
        
        filename = regexprep(plot_type, '\s+', '');
        
        fig = figure('numbertitle','off','name',plot_type,'visible','off');
        if(log_flag)
            halllines = semilogy(t,all_traj,'Color',seccolor(1,:));
        else
            halllines = plot(t,all_traj,'Color',seccolor(1,:));
        end
        hallgroup = hggroup;
        set(halllines,'Parent',hallgroup);
        set(get(get(hallgroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
        hold on
        if(log_flag)
            semilogy(t,mean_traj,'Color',maincolor(1,:)*.7,'LineWidth',2)
            semilogy(t,ODE_traj,'--','Color',maincolor(1,:),'LineWidth',2)
        else
            plot(t,mean_traj,'Color',maincolor(1,:)*.7,'LineWidth',2)
            plot(t,ODE_traj,'--','Color',maincolor(1,:),'LineWidth',2)
        end
        xlabel('Weeks', 'fontsize',fsize)
        ylabel(plot_type, 'fontsize',fsize)
        set(gca,'fontsize',fsize);
        legend('All Gillespie Runs','Mean of Runs','ODE Run')
        print(fig,'-djpeg',[filepath filename '.jpg']);
        close(fig);
        
        fig = figure('numbertitle','off','name',[plot_type ' PDF'],'visible','off');
        count = 1;
        Xdata = zeros(t_steps,100);
        Ydata = zeros(t_steps,100);
        for j=1:5:(length(t)-9)
            temp = reshape(all_traj(j:j+9,:),[],1);
            [Ydata(count,:), Xdata(count,:)] = hist(temp,100);
            Ydata(count,:) = (Ydata(count,:)/max(Ydata(count,:)))*10 + t(j);
            count = count + 1;
        end
        if(log_flag)
            h_PDFlines = semilogy(Ydata',Xdata','color',maincolor(1,:));
        else
            h_PDFlines = plot(Ydata',Xdata','color',maincolor(1,:));
        end
        h_PDFgroup = hggroup;
        set(h_PDFlines,'Parent',h_PDFgroup);
        set(get(get(h_PDFgroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
        hold on
        if(log_flag)
            semilogy(t,mean_traj,'Color',maincolor(1,:)*.7,'LineWidth',2)
            semilogy(t,ODE_traj,'--','Color',maincolor(1,:),'LineWidth',2)
        else
            plot(t,mean_traj,'Color',maincolor(1,:)*.7,'LineWidth',2)
            plot(t,ODE_traj,'--','Color',maincolor(1,:),'LineWidth',2)
        end
        xlabel('Weeks', 'fontsize',fsize)
        ylabel(plot_type, 'fontsize',fsize)
        set(gca,'fontsize',fsize);
        legend('PDF of Gillespie','Mean of Runs','ODE Run')
        saveas(fig,[filepath filename '_PDF'],'fig');
        print(gcf,'-djpeg',[filepath filename '_PDF.jpg']);
        close(fig);
        
        fig = figure('numbertitle','off','name',[plot_type ' w/error bars'],'visible','off');
        if(log_flag)
            semilogy(t,mean_traj,'Color',maincolor(1,:)*.7,'LineWidth',2)
        else
            plot(t,mean_traj,'Color',maincolor(1,:)*.7,'LineWidth',2)
        end
        hold on
        if(log_flag)
            semilogy(t,ODE_traj,'--','Color',maincolor(1,:),'LineWidth',2)
        else
            plot(t,ODE_traj,'--','Color',maincolor(1,:),'LineWidth',2)
        end
        legend('Mean of Runs','ODE Run')
        %warning('off','MATLAB:Axes:NegativeDataInLogAxis');
        if(log_flag)
            errorbar(t,mean_traj,mean_traj-exp(log(mean_traj)-std(log(all_traj),0,2)),mean_traj-exp(log(mean_traj)+std(log(all_traj),0,2)),'kx')
        else
            errorbar(t,mean_traj,std(all_traj,0,2),'kx')
        end
        xlabel('Weeks', 'fontsize',fsize)
        ylabel(plot_type, 'fontsize',fsize)
        xlim([min(t), max(t)]);
        set(gca,'fontsize',fsize);
        saveas(fig,[filepath filename '_EB'],'fig');
        print(gcf,'-djpeg',[filepath filename '_EB.jpg']);
        %warning('on','MATLAB:Axes:NegativeDataInLogAxis');
        close(fig);
        
    end

end

%% Calc P1 and P2 trajectories gilsp

function [P1_traj, P2_traj] = P_calc(n_snap, n_P1, n_P2, traj_array, resIdx, nRes, start, finish, n_pop, P2traj)
P1_traj = zeros(n_snap,n_P1);
P2_traj = zeros(n_snap,n_P2);

uniqueResIdx = 1:21;              % unique residue indices residing in MSA;

parfor snap=1:n_snap
    P1_temp = zeros(1,n_P1);
    P2_temp = zeros(1,n_P2);
    
    % one-site residue ordering and target probabilities
    hist1=histc(squeeze(traj_array(snap,:,:)),uniqueResIdx);              % one-site mutational histograms
    %    [hist1,~]=sort(hist1,1,'descend');     % reordering the residue indices at each site in descending order of prevalence
    
    P1_model=cell(1,(finish-start+1));
    for i=1:(finish-start+1)
        P1_model{i}=hist1(resIdx{start+i-1},i)/n_pop; %#ok<PFBNS>
    end
    
    
    % two-site target probabilities
    % \-> matrix indexing based on the one-site residue ordering in each pair
    if(P2traj || snap == n_snap)
        P2_model=P2_calc(traj_array(snap,:,:),resIdx(start:finish),(finish-start+1),n_pop);
        
        pos = 1;
        for i=1:(finish-start+1)
            for j=i+1:(finish-start+1)
                for q=1:nRes(start+j-1)  %#ok<PFBNS>
                    P2_temp(1,pos:(pos+(nRes(start+i-1))-1)) = P2_model{i,j}(:,q);
                    pos = pos + nRes(start+i-1);
                end
            end
        end
        P2_traj(snap,:) = P2_temp;
        
        
    end
    
    pos = 1;
    for i=1:(finish-start+1)
        P1_temp(1,pos:(pos+nRes(start+i-1))-1) = P1_model{i};
        pos = pos + nRes((start+i-1));
    end
    P1_traj(snap,:) = P1_temp;
    
end
end

%% Calculate P2
% This function is needed to allow P_calc to use parfor

function P2_model=P2_calc(traj_array,resIdx,m,n_pop)
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
                p=hist3(squeeze(traj_array(1,:,[i,j])),'Edges',{edge_i,edge_j});
            else % length(edge_j)==1
                p=hist3(squeeze(traj_array(1,:,[i,j])),'Edges',{edge_i,cat(1,edge_j,200184*edge_j)});
                p=p(:,1);
            end
        else	% length(edge_i)==1
            if ( length(edge_j)>1 )
                p=hist3(squeeze(traj_array(1,:,[i,j])),'Edges',{cat(1,edge_i,170786*edge_i),edge_j});
                p=p(1,:);
            else % length(edge_j)==1
                p=hist3(squeeze(traj_array(1,:,[i,j])),'Edges',{cat(1,edge_i,170786*edge_i),cat(1,edge_j,200184*edge_j)});
                p=p(1,1);
            end
        end
        
        p=p(idx_unsort_i,idx_unsort_j)/n_pop;                     % reordering and normalizing
        P2_model{i,j}=p;
    end
end
end
