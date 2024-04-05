% load and plot iemmod results
% This script is intended for the distribution in the open-source package
clear all; close all; clc; tic; 

doplotmisfit = 1; % Plot misfit (1) or not (0)
doplotmod = -1; % < 0: Plot the last model
                % = 0: Do not plot any model
                % > 0: Plot the model of the specified iteration number
doprint = 0; % Print figures (1) or no (0)
                
z = [0,200,1200,1220]; % depth of layer interfaces
folder = '/data_a/hunziker/testing/rTE_only/'; % folder where inversion results are stored
name = 'R^{TE}';
truemod = [1,1,1,0.02,1]; % conductivities of the true model
startmod = [1,1,0.5,0.5,0.5]; % conductivities of the starting model
filestem = 'econd'; % stem of the filename containing the inversion results
startz = -100; % the plot starts at this depth
endz = 1500; % the plot ends at this depth
fs = 16; % Fontsize

% Determine how many iterations are completed
checker = 1;
in = 1;
while checker ~= 0
    checker = exist([deblank(folder),filestem,'_',num2str(in),'.bin'],'file');
    in = in + 1;
end
nof = in-2;
fprintf('Completed iterations: %d\n',nof);

% Determine how many subsurface parameters are in the model
fid = fopen([deblank(folder),filestem,'_1.bin'],'r');
status = fseek(fid,0,'eof');
filesize = ftell(fid);
nel = round(filesize/8); % this contains all the conductivities plus the misfit
fclose(fid);
fprintf('Amount of subsurface parameters: %d\n',nel-1);

% Load the model and the misfit of each iteration
data = zeros(nof,nel-1);
misfit = zeros(nof,1);
for in = 1:nof
    fid = fopen([deblank(folder),filestem,'_',num2str(in),'.bin'],'r');
    temp = fread(fid,nel,'float64',0,'ieee-le');
    fclose(fid);
    data(in,1:nel-1) = temp(1:nel-1);
    misfit(in) = temp(nel);
end

% Plot the misfit as a function of iterations
if doplotmisfit==1
    figure(1);
    semilogy(misfit,'k');
    grid off
    if doplotmod>0
        % Plot the convergence history only until the selected iteration
        xlim([1 doplotmod]);
    else
        % Plot the complete available convergence history
        xlim([1 nof])
    end
    xlabel('Iteration','Fontsize',fs);
    ylabel('Misfit','Fontsize',fs);
    title(['Convergence ',deblank(name)],'Fontsize',fs);
    set(gca,'Fontsize',fs);
    if doprint == 1
        print('-dtiff','convergence.tif');
    end
end

if doplotmod>nof
    error(['Wanted to plot results after iteration ',num2str(doplotmod),', but found only ',num2str(nof),' Iterations.']);
end
if doplotmod~=0
    % Set up the depth vector for plotting
    zvecplot = zeros(2*length(z)+2,1);
    zvecplot(1) = startz;
    iz = 1;
    for in = 2:2:length(zvecplot)-1
        zvecplot(in) = z(iz);
        zvecplot(in+1) = z(iz);
        iz = iz+1;
    end
    zvecplot(end) = endz;
    
    % Set up the conductivity vector for plotting
    startmodplot = zeros(2*length(startmod),1);
    truemodplot = zeros(2*length(startmod),1);
    selmodplot = zeros(2*length(startmod),1);
    if doplotmod<0
        % The selected model for plotting is the latest available model
        selmod = nof;
    else
        % The selected model for plotting is specified in the variable doplotmod
        selmod = doplotmod;
    end
    iz = 1;
    for in = 1:2:length(startmodplot)
        startmodplot(in) = startmod(iz);
        startmodplot(in+1) = startmod(iz);
        truemodplot(in) = truemod(iz);
        truemodplot(in+1) = truemod(iz);
        selmodplot(in) = data(selmod,iz);
        selmodplot(in+1) = data(selmod,iz);
        iz = iz+1;
    end
    
    % Plot the model
    figure(2);
    plot(truemodplot,zvecplot,'b','Linewidth',2);
    hold on;
    plot(startmodplot,zvecplot,'g');
    plot(selmodplot,zvecplot,'r');
    hold off;
    xlim([-0.1 3.1])
    ylim([startz endz])
    xlabel('Conductivity [S/m]','Fontsize',fs);
    ylabel('Depth [m]','Fontsize',fs);
	legend('True Model','Starting Model','Estimated Model','Location','SouthEast');
    title(['Model ',deblank(name)],'Fontsize',fs);
    set(gca,'YDir','reverse','Fontsize',fs);
    if doprint == 1
        print('-dtiff','models.tif');
    end
end

toc