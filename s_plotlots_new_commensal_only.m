%%%%%Code used to produce Fig 2 

function s_plotlots_new_commensal_only
% Read in results and plot various things that might be useful :)

% Set parameters
tsph=10;        % time steps per hour
tth=1440;         % total time in hours
tt=tth*tsph;
plothours=[0 24 48 72 96 720];
% plothours=[0 12 24 48 72 96 720];

numplots=6;     % Number of plots to make across time period
Ks = 1e9;     % Carrying capacity in small intestine
Kl = 1e12;    % Carrying capacity in large intestine
n_pat = 60;     % number of patches

% Set up filename
fileA='Bacteria_n_only_Ks_';
fileB=sprintf('%d',Ks);
fileBX='_Kl_';
fileBY=sprintf('%d',Kl);
fileC='_patches_';
fileD=sprintf('%d',n_pat);
fileE='_normal';
fileF='_new';
fileE1='_inflammation-causing';
fileE2='_inflammation-level';
fileX='_end';
fileY=sprintf('%d',tt);
fileG='.mat';

fileM1={fileA,fileB,fileBX,fileBY,fileC,fileD,fileX,fileY,fileE,fileG};
filename1=strjoin(fileM1,'');

load(filename1)

fileM2={fileA,fileB,fileBX,fileBY,fileC,fileD,fileX,fileY,fileE1,fileG};
filename2=strjoin(fileM2,'');

load(filename2)

fileM3={fileA,fileB,fileBX,fileBY,fileC,fileD,fileX,fileY,fileE2,fileG};
filename3=strjoin(fileM3,'');

load(filename3)

fileM4={fileA,fileB,fileBX,fileBY,fileC,fileD,fileX,fileY,fileF,fileG};
filename4=strjoin(fileM4,'');

load(filename4)

% Set up plots
%plothours=linspace(1,tt+1,numplots);
plothours=plothours*tsph+1;

spaceplot=linspace(1,n_pat,n_pat);
% 


figure(1);
clf;
subplot(2,3,1)
for i=1:numplots
    time=plothours(i);
    ynorm=x_normal(time,:);
    semilogy(spaceplot,ynorm,'Linewidth',2)
    title('Old Commensal Bacteria')
    xlabel('Patches'); ylabel('Bacterial load')
     legend('t=0h', 't=12h','t=48h', 't=72h', 't=96h', 't=720h')
    set(gca, 'Linewidth',2,'Fontsize',18)
    hold on
end


subplot(2,3,2)
for i=1:numplots
    time=plothours(i);
    ynew=x_new(time,:);
    semilogy(spaceplot,ynew,'Linewidth',2)
    title(' New Commensal Bacteria')
    xlabel('Patches');  ylabel('Bacterial load')
    legend('t=0h', 't=6h', 't=12h', 't=24h', 't=48h', 't=72h', 't=96h', 't=720h')
    set(gca, 'Linewidth',2,'Fontsize',18)
    hold on
    
end


subplot(2,3,3)
for i=1:numplots
    time=plothours(i);
    ycommensal=x_new(time,:)+x_normal(time,:);
    semilogy(spaceplot,ycommensal,'Linewidth',2)
    title('All Commensal Bacteria')
    xlabel('Patches'); ylabel('Bacterial load')
    legend('t=0h', 't=6h', 't=12h', 't=24h', 't=48h', 't=72h', 't=96h', 't=720h')
    set(gca, 'Linewidth',2,'Fontsize',18)
    hold on
    
end
% plothours=[0 48 120 360 720 1440];
% plothours=plothours*tsph+1;
% time=plothours(2);
% ynorm=log10(x_normal(time,:));
% ynorm(isinf(ynorm)) = 0; 
% ynew=log10(x_new(time,:));
% ynew(isinf(ynew)) = 0; 
% 
% time1=plothours(4);
% ynorm1=log10(x_normal(time1,:));
% ynorm1(isinf(ynorm1)) = 0; 
% yres1=log10(x_new(time1,:));
% yres1(isinf(yres1)) = 0; 
% 
% time2=plothours(5);
% ynorm2=log10(x_normal(time2,:));
% ynorm2(isinf(ynorm2)) = 0; 
% ynew2=log10(x_new(time2,:));
% ynew2(isinf(ynew2)) = 0; 
% 
% time3=plothours(6);
% ynorm3=log10(x_normal(time3,:));
% ynorm3(isinf(ynorm3)) = 0; 
% ynew3=log10(x_new(time3,:));
% ynew3(isinf(ynew3)) = 0; 
% 
% X=[ynorm ynew]';X1=[ynorm1 yres1]';X2=[ynorm2 ynew2]';X3=[ynorm3 ynew3]';
% 
% Y1=repmat('Old',length(ynorm),1);
% Y2=repmat('New',length(ynorm),1);
% Bacteria=[Y1;Y2];
% Bacteria = cellstr(Bacteria);
% 
% grouporder={'Old','New'};
% subplot(2,3,4)
% vs = violinplot(X, Bacteria,'GroupOrder',grouporder,'Width',0.2,'ViolinAlpha',0.4);
% ylabel('Log_{10} Bacterial density');
% title('Old and new commensal at day 2')
% set(gca, 'Linewidth',2,'Fontsize',18)
% 
% subplot(2,3,5)
% vs1 = violinplot(X1,Bacteria,'GroupOrder',grouporder,'Width',0.2,'ViolinAlpha',0.4);
% ylabel('Log_{10} Bacterial density');
% title('Old and new commensal at day 15')
% set(gca, 'Linewidth',2,'Fontsize',18)





% 
% 
% %%%%%%%Contour plots
% figure(4);
% clf;
% subplot(2,2,1)
% pcolor(spaceplot,t/240,log10(x_normal));
% shading interp;hold on;
% [M,c]=contour(spaceplot,t/240,log10(x_normal),10);
% c.LineColor = 'k';
% colormap parula
% h=colorbar;
% set(get(h,'label'),'string','Log_{10} Bacterial density');
% xlabel('Patch'); ylabel('Time (days)');
% % set(gca, 'clim', [min(Ci(:)) max(Ci(:))])
% title('Old Commensal Bacteria')
% set(gca, 'Linewidth',2,'Fontsize',18)
% 
% subplot(2,2,2)
% pcolor(spaceplot,t/240,log10(x_new));
% shading interp;hold on;
% [M,c]=contour(spaceplot,t/240,log10(x_new),10);
% c.LineColor = 'k';
% colormap parula
% h=colorbar;
% set(get(h,'label'),'string','Log_{10} Bacterial density');
% xlabel('Patch'); ylabel('Time (days)');
% title('New Commensal Bacteria')
% set(gca, 'Linewidth',2,'Fontsize',18)
% 
% subplot(2,2,3)
% pcolor(spaceplot,t/240,log10(x_normal+x_new));
% shading interp;hold on;
% [M,c]=contour(spaceplot,t/240,log10(x_normal+x_new),10);
% c.LineColor = 'k';
% colormap parula
% h=colorbar;
% set(get(h,'label'),'string','Log_{10} Bacterial density');
% xlabel('Patch'); ylabel('Time (days)');
% % set(gca, 'clim', [min(Ci(:)) max(Ci(:))])
% title('All Commensal Bacteria')
% set(gca, 'Linewidth',2,'Fontsize',18)








plothours=[0 48 120 360 720 1440];
plothours=plothours*tsph+1;
time=plothours(2);
ynorm=log10(x_normal(time,:));
ynorm(isinf(ynorm)) = 0; 
ynew=log10(x_new(time,:));
ynew(isinf(ynew)) = 0; 

time1=plothours(4);
ynorm1=log10(x_normal(time1,:));
ynorm1(isinf(ynorm1)) = 0; 
yres1=log10(x_new(time1,:));
yres1(isinf(yres1)) = 0; 

time2=plothours(5);
ynorm2=log10(x_normal(time2,:));
ynorm2(isinf(ynorm2)) = 0; 
ynew2=log10(x_new(time2,:));
ynew2(isinf(ynew2)) = 0; 

time3=plothours(6);
ynorm3=log10(x_normal(time3,:));
ynorm3(isinf(ynorm3)) = 0; 
ynew3=log10(x_new(time3,:));
ynew3(isinf(ynew3)) = 0; 

X=[ynorm ynew]';X1=[ynorm1 yres1]';X2=[ynorm2 ynew2]';X3=[ynorm3 ynew3]';

Y1=repmat('Old',length(ynorm),1);
Y2=repmat('New',length(ynorm),1);
Bacteria=[Y1;Y2];
Bacteria = cellstr(Bacteria);

figure(11);
clf;

subplot(2,2,1)
vs = violinplot(X, Bacteria,'Width',0.2,'ViolinAlpha',0.4);
ylabel('Log_{10} Bacterial density');
title('Old and new commensal at day 2')

subplot(2,2,2)
vs1 = violinplot(X1,Bacteria,'Width',0.2,'ViolinAlpha',0.4);
ylabel('Log_{10} Bacterial density');
title('Old and new commensal at day 15')

subplot(2,2,3)
vs2 = violinplot(X2,Bacteria,'Width',0.2,'ViolinAlpha',0.4);
ylabel('Log_{10} Bacterial density');
title('Old and new commensal at day 30')

subplot(2,2,4)
vs3 = violinplot(X3, Bacteria,'Width',0.2,'ViolinAlpha',0.4 );
ylabel('Log_{10} Bacterial density');
title('Old and new commensal at day 60')


end

