%%%%%%%%%Code used to generate Fig S7

function s_plotlots_appendiceal
% Set parameters
tsph=10;        % time steps per hour
tth=720;         % total time in hours
tt=tth*tsph;
plothours=[0 12 48 72 96 720];
% plothours=[0 12 24 48 72 96 720];

numplots=6;     % Number of plots to make across time period
Ks = 1e9;     % Carrying capacity in small intestine
Kl = 1e12;    % Carrying capacity in large intestine
n_pat = 60;     % number of patches
Ka = 1e6;  
% Set up filename
fileA='Bacteria_Ks_';
% fileA='Bacteria_n_only_Ks_';
fileB=sprintf('%d',Ks);
fileBX='_Kl_';
fileBY=sprintf('%d',Kl);
fileBXa='_Ka_';
fileBYa=sprintf('%d',Ka);
fileC='_patches_';
fileD=sprintf('%d',n_pat);
fileE='_normal';
fileF='_resistant';
fileE1='_inflammation-causing';
fileE2='_inflammation-level';
fileE3='_appendiceal';
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


fileM5={fileA,fileBXa,fileBYa,fileX,fileY,fileE3,fileG};
filename5=strjoin(fileM5,'');

load(filename5)


% Set up plots
%plothours=linspace(1,tt+1,numplots);
% 


%%%%%%%Contour plots

spaceplot=linspace(1,n_pat,n_pat);
figure(4);
clf;
subplot(1,2,1)
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
title('Commensal Bacteria')
set(gca, 'Linewidth',2,'Fontsize',18)



subplot(1,2,2)
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

