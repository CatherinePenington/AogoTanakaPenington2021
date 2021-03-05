%%%%%Code used to produce Fig S2 
function s_plotlots_pathogen_day
% Read in results and plot various things that might be useful :)

% Set parameters
tsph=10;        % time steps per hour
tth=720;         % total time in hours
tt=tth*tsph;
% plothours=[0 6 12 24 48 72 96 720];
plothours=[0 72 120 240 360 480 600 720];

numplots=8;     % Number of plots to make across time period
Ks = 1e9;     % Carrying capacity in small intestine
Kl = 1e12;    % Carrying capacity in large intestine
n_pat = 60;     % number of patches

% Set up filename
fileA='Bacteria_Ks_';
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
subplot(2,2,1)
for i=1:numplots
    time=plothours(i);
    ynorm=x_normal(time,:);
    semilogy(spaceplot,ynorm,'Linewidth',2)
    title('Old Commensal Bacteria')
    xlabel('Patch'); ylabel('Bacterial density')
%      legend('t=0h', 't=6h', 't=12h', 't=24h', 't=48h', 't=72h', 't=96h', 't=720h')
     legend('t=0h', 't=3d', 't=5d', 't=10d', 't=15d', 't=20d', 't=25d', 't=30d')
    set(gca, 'Linewidth',2,'Fontsize',18)
    hold on
end


subplot(2,2,2)
for i=1:numplots
    time=plothours(i);
    ynew=x_new(time,:);
    semilogy(spaceplot,ynew,'Linewidth',2)
    title(' New Commensal Bacteria')
    xlabel('Patch');  ylabel('Bacterial density')
%     legend('t=0h', 't=6h', 't=12h', 't=24h', 't=48h', 't=72h', 't=96h', 't=720h')
    legend('t=0h', 't=3d', 't=5d', 't=10d', 't=15d', 't=20d', 't=25d', 't=30d')
    set(gca, 'Linewidth',2,'Fontsize',18)
    hold on
    
end


subplot(2,2,3)
for i=1:numplots
    time=plothours(i);
    ycommensal=x_new(time,:)+x_normal(time,:);
    semilogy(spaceplot,ycommensal,'Linewidth',2)
    title('All Commensal Bacteria')
    xlabel('Patch'); ylabel('Bacterial density')
%     legend('t=0h', 't=6h', 't=12h', 't=24h', 't=48h', 't=72h', 't=96h', 't=720h')
    legend('t=0h', 't=3d', 't=5d', 't=10d', 't=15d', 't=20d', 't=25d', 't=30d')
    set(gca, 'Linewidth',2,'Fontsize',18)
    hold on
    
end

subplot(2,2,4)

for i=1:numplots
    time=plothours(i);
    yinfl=x_infl(time,:);
    semilogy(spaceplot,yinfl,'Linewidth',2)
    title('Inflammatory Bacteria')
    xlabel('Patch'); ylabel('Bacterial density')
%      legend('t=0h', 't=6h', 't=12h', 't=24h', 't=48h', 't=72h', 't=96h', 't=720h')
     legend('t=0h', 't=3d', 't=5d', 't=10d', 't=15d', 't=20d', 't=25d', 't=30d')
    set(gca, 'Linewidth',2,'Fontsize',18)
    hold on
end

%%%inflammatory response
figure(2);
clf;
for i=1:numplots
    time=plothours(i);
    yil=Il(time,:);
    plot(spaceplot,yil,'Linewidth',2)
    title('Inflammatory response')
    xlabel('Patch'); ylabel('Inflammation level')
%     legend('t=0h', 't=6h', 't=12h', 't=24h', 't=48h', 't=72h', 't=96h', 't=720h')
    legend('t=0h', 't=3d', 't=5d', 't=10d', 't=15d', 't=20d', 't=25d', 't=30d')
    set(gca, 'Linewidth',2,'Fontsize',18)
    hold on
    
end



%%%%%%%Contour plots
figure(4);
clf;
subplot(2,2,1)
pcolor(spaceplot,t/240,log10(x_normal));
shading interp;hold on;
[M,c]=contour(spaceplot,t/240,log10(x_normal),10);
c.LineColor = 'k';
colormap parula
h=colorbar;
caxis([0 12]);
set(get(h,'label'),'string','Log_{10} Bacterial density');
xlabel('Patch'); ylabel('Time (days)');
% set(gca, 'clim', [min(Ci(:)) max(Ci(:))])
title('Old Commensal Bacteria')
set(gca, 'Linewidth',2,'Fontsize',18)

subplot(2,2,2)
pcolor(spaceplot,t/240,log10(x_new));
shading interp;hold on;
[M,c]=contour(spaceplot,t/240,log10(x_new),10);
c.LineColor = 'k';
colormap parula
h=colorbar;
caxis([0 12]);
set(get(h,'label'),'string','Log_{10} Bacterial density');
xlabel('Patch'); ylabel('Time (days)');
title('New Commensal Bacteria')
set(gca, 'Linewidth',2,'Fontsize',18)

subplot(2,2,3)
pcolor(spaceplot,t/240,log10(x_normal+x_new));
shading interp;hold on;
[M,c]=contour(spaceplot,t/240,log10(x_normal+x_new),10);
c.LineColor = 'k';
colormap parula
h=colorbar;
caxis([0 12]);
set(get(h,'label'),'string','Log_{10} Bacterial density');
xlabel('Patch'); ylabel('Time (days)');
% set(gca, 'clim', [min(Ci(:)) max(Ci(:))])
title('All Commensal Bacteria')
set(gca, 'Linewidth',2,'Fontsize',18)

subplot(2,2,4)
pcolor(spaceplot,t/240,log10(x_infl));
shading interp;hold on;
[M,c]=contour(spaceplot,t/240,log10(x_infl),10);
c.LineColor = 'k';
colormap parula
h=colorbar;
caxis([0 12]);
set(get(h,'label'),'string','Log_{10} Bacterial density');
xlabel('Patch'); ylabel('Time (days)');
title('Inflammatory Bacteria')
set(gca, 'Linewidth',2,'Fontsize',18)


%%%%%%%Contour plot Inflammation
figure(5);
pcolor(spaceplot,t/240,Il);
shading interp;hold on;
[M,c]=contour(spaceplot,t/240,Il,10);
c.LineColor = 'k';
colormap parula
h=colorbar;
caxis([0 0.7]);
set(get(h,'label'),'string','Inflammation level');
xlabel('Patch'); ylabel('Time (days)');
title('Inflammatory response')
set(gca, 'Linewidth',2,'Fontsize',18)


end

