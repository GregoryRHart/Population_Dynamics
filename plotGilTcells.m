function plotGilTcells(filepath)
% This function reads in data from Tcell_traj.dat from multiple runs of the
% Gillespie version of the Population Dynamics code and from a run of the
% ODE Population Dynamics code. The function goes on to make plots of this
% data to compare the ODE and mean behavior as well as get a feeling for
% the varience. For more information on the Tcell_traj.dat file see
% load_TcellDat.m
%
% VAR IN
% filepath  - filepath is an optional input arguement that allows the user
%             to specifically the parent directory of all the individual
%             runs. If it is not given the function reads from the current
%             working directory.
%
% VAR OUT
% none
%
% FILE IN
% Tcell_traj.dat - This script reads in Tcell_traj.dat from subdirectories
%                  that are numeric ([0-9]*) and from ODE. Tcell_traj.dat
%                  is output by the Population Dynamics code. The numeric
%                  directories should contain runs of the Gillespie version
%                  of the code with different random seeds (preferably the
%                  same as the name of the directory) and the ODE
%                  directory should contain a run of the ODE version of the
%                  code. Tcell_traj.dat will be empty if the code was run
%                  with no immune pressure.
%
% FILE OUT
% Naive_cells           - Plot of the respective cell type versus time. The
% Effector_cells          trajectories from all of the individual runs are
% Memory_cells            plotted in washed thin lines. Then the mean
% Susceptible_cells       trajectory and the ODE trajectory are ploted in
%                         thicker darker lines. Only a .jpg is saved.
%
% Naive_cells_PDF       - Plot of the respective cell type versus time. The
% Effector_cells_PDF      mean trajectory and the ODE trajectory are plotted
% Memory_cells_PDF        and then the data from individual runs are
% Susceptible_cells_PDF   aggrogated over ten time step windows to generate
%                         PDFs of the trajectory which is plotted on top of
%                         the existing graph. Both a .fig and .jpg are
%                         saved.
%
% Naive_cells_ED        - Plot of the respective cell type versus time. The
% Effector_cells_ED       mean trajectory and the ODE trajectory are plotted
% Memory_cells_ED         and then the individual runs are used to
% Susceptible_cells_ED    calculate error bars which are added to the plot.
%                         Both a .fig and .jpg are saved.
%
% Written by:
% GRH 21 Feb 2017
%

%% Setup

tic;

if nargin < 1
    filepath ='.';
end
if(~isempty(filepath) && filepath(end) ~= filesep)
    filepath=[filepath filesep];
end

fsize=12;

d = dir(filepath);
d = {d([d.isdir] == 1).name};
d = regexpi(d,'[0-9]*','match');
d = [d{:}];

%% Read in data

% get size
[~, n_epitope, n_snap, ~, traj_array] = evalc('load_TcellDat([filepath d{1}])');

% preallocate matrices
Nall = zeros(n_snap,n_epitope,length(d));
Eall = zeros(n_snap,n_epitope,length(d));
Mall = zeros(n_snap,n_epitope,length(d));
Iall = zeros(n_snap,n_epitope,length(d));

% Add data from the first directory
Nall(:,:,1) = traj_array(:,1:4:end);
Eall(:,:,1) = traj_array(:,2:4:end);
Mall(:,:,1) = traj_array(:,3:4:end);
Iall(:,:,1) = traj_array(:,4:4:end);

% read in data for from the remaining runs
for k=2:length(d)
    [~, ~, ~, ~, traj_array] = evalc('load_TcellDat([filepath d{k}])');
    Nall(:,:,k) = traj_array(:,1:4:end);
    Eall(:,:,k) = traj_array(:,2:4:end);
    Mall(:,:,k) = traj_array(:,3:4:end);
    Iall(:,:,k) = traj_array(:,4:4:end);
end

% find mean behavoir
N = mean(Nall,3);
E = mean(Eall,3);
M = mean(Mall,3);
I = mean(Iall,3);

% read in ODE solution
[~,~, ~, t, traj_array] = evalc('load_TcellDat([filepath ''ODE/''])');

NODE = traj_array(:,1:4:end);
EODE = traj_array(:,2:4:end);
MODE = traj_array(:,3:4:end);
IODE = traj_array(:,4:4:end);

%% Plotting

maincolor = [0,0,1; 0,1,0; 1,0,0; 0,1,1; 1,0,1; 1,1,0; 0,0,0];
seccolor = maincolor*.5+.5;
t_steps = ceil((length(t)-9)/5);
% Get name of epitopes
epitopes = importdata([filepath d{1} '/epitopes.dat']);
[~, epitopes] = strtok(epitopes, '_');
[~, epitopes] = strtok(epitopes, '_');
[epitopes,~] = strtok(epitopes,'.');
for k=1:n_epitope
    epitopes{k} = epitopes{k}(2:end);
    epitopes{k} = strrep(epitopes{k},'_',' ');
end

plot_cellTraj(Nall, N, NODE, 'Naive');
plot_cellTraj(Eall, E, EODE, 'Effector');
plot_cellTraj(Mall, M, MODE, 'Memory');
plot_cellTraj(Iall, I, IODE, 'Susceptible');

toc;

    % This function makes three different plots of the trajectory for a
    % specific type of T-cell.
    function plot_cellTraj(all_traj, mean_traj, ODE_traj, cell_type)
        % Plot the trajectory from all runs
        hTcelllines = cell(n_epitope,1);
        figure('numbertitle','off','name',[cell_type ' Cells'],'visible','off')
        for i=1:n_epitope
            hTcelllines{i} = plot(t,squeeze(all_traj(:,i,:)),'Color',seccolor(i,:));
            hold on
        end
        hTcelllines = [hTcelllines{:}];
        hTcellgroup = hggroup;
        set(hTcelllines,'Parent',hTcellgroup);
        set(get(get(hTcellgroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
        for i=1:n_epitope
            plot(t,mean_traj(:,i),'Color',maincolor(i,:)*.7,'LineWidth',2)
            plot(t,ODE_traj(:,i),'--','Color',maincolor(i,:),'LineWidth',2)
        end
        xlabel('Weeks', 'fontsize',fsize)
        ylabel('# of cells', 'fontsize',fsize)
        set(gca,'fontsize',fsize);
        label = cell(n_epitope*2+1,1);
        label{1} = 'All runs';
        for i=1:n_epitope
            label{(i-1)*2+2} = [epitopes{i}, ' mean'];
            label{(i-1)*2+3} = [epitopes{i}, ' ODE'];
        end
        legend(label{:})
        print(gcf,'-djpeg',[filepath cell_type '_cells.jpg']);
        close(gcf);
        
        % Plot the mean trajectory and ODE trajectory with the PDF of the
        % runs.
        y = [min(ODE_traj(1,:)),max(ODE_traj(1,:))];
        fig = figure('numbertitle','off','name',[cell_type ' Cells PDF'],'visible','off');
        for i=1:n_epitope
            Trange = [floor(min(min(min(all_traj(:,i,:))))),ceil(max(max(max(all_traj(:,i,:)))))];
            nbins = min(100, Trange(2)-Trange(1) + 3);
            count = 1;
            Xdata = zeros(t_steps,nbins);
            Ydata = zeros(t_steps,nbins);
            for j=1:5:(length(t)-9)
                temp = reshape(all_traj(j:j+9,i,:),[],1);
                [Ydata(count,:), Xdata(count,:)] = hist(temp,nbins);
                Ydata(count,:) = (Ydata(count,:)/max(Ydata(count,:)))*10 + t(j);
                count = count + 1;
            end
            hTcelllines = plot(Ydata',Xdata','color',maincolor(i,:));
            hold on
            y = [min([y,min(Xdata)]),max([y,max(Xdata)])];
            hTcellgroup = hggroup;
            set(hTcelllines,'Parent',hTcellgroup);
            set(get(get(hTcellgroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','on');
            plot(t,mean_traj(:,i),'Color',maincolor(i,:)*.7,'LineWidth',2)
            plot(t,ODE_traj(:,i),'--','Color',maincolor(i,:),'LineWidth',2)
        end
        
        y = [min([y, min(mean_traj), min(ODE_traj)]), max([y, max(mean_traj), max(ODE_traj)])];
        ylim(y)
        
        xlabel('Weeks', 'fontsize',fsize)
        ylabel('# of cells', 'fontsize',fsize)
        set(gca,'fontsize',fsize);
        label = cell(n_epitope*3,1);
        label{1} = 'All runs';
        for i=1:n_epitope
            label{(i-1)*3+1} = [epitopes{i}, ' PDF'];
            label{(i-1)*3+2} = [epitopes{i}, ' mean'];
            label{(i-1)*3+3} = [epitopes{i}, ' ODE'];
        end
        legend(label{:})
        saveas(fig,[filepath cell_type '_cells_PDF'],'fig');
        print(fig,'-djpeg',[filepath cell_type '_cells_PDF.jpg']);
        close(fig);
        
        % Plot the mean trajectory and ODE trajectory with error bars.
        fig = figure('numbertitle','off','name',[cell_type ' Cells w/error bars'],'visible','off');
        plot(t,mean_traj(:,i),'Color',maincolor(1,:)*.7,'LineWidth',2)
        hold on
        plot(t,ODE_traj(:,i),'--','Color',maincolor(1,:),'LineWidth',2)
        for i=2:n_epitope
            plot(t,mean_traj(:,i),'Color',maincolor(i,:)*.7,'LineWidth',2)
            hold on
            plot(t,ODE_traj(:,i),'--','Color',maincolor(i,:),'LineWidth',2)
        end
        label = cell(n_epitope*2,1);
        for i=1:n_epitope
            label{(i-1)*2+1} = [epitopes{i}, ' mean'];
            label{(i-1)*2+2} = [epitopes{i}, ' ODE'];
        end
        legend(label{:})
        for i=1:n_epitope
            errorbar(t,mean_traj(:,i),std(squeeze(all_traj(:,i,:)),0,2),'x','color',maincolor(i,:))
        end
        xlabel('Weeks', 'fontsize',fsize)
        ylabel('# of cells', 'fontsize',fsize)
        xlim([min(t), max(t)]);
        set(gca,'fontsize',fsize);
        saveas(fig,[filepath cell_type '_cells_EB'],'fig');
        print(gcf,'-djpeg',[filepath cell_type '_cells_EB.jpg']);
        close(fig)
    end

end
