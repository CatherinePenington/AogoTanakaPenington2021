%%%%Code used to generate Fig 7 & 8
function s_plotlots_resistant_initially
% Set parameters
tsph=10;        % time steps per hour
tth=1440;         % total time in hours
tt=tth*tsph;
% plothours=[0 6 12 24 48 72 96 720];
plothours=[0 72 120 240 360 480 600 720];

numplots=8;     % Number of plots to make across time period
Ks = 1e9;     % Carrying capacity in small intestine
Kl = 1e12;    % Carrying capacity in large intestine
n_pat = 60;     % number of patches

% Set up filename
fileA='Bacteria_Ks_';
% fileA='Bacteria_n_only_Ks_';
fileB=sprintf('%d',Ks);
fileBX='_Kl_';
fileBY=sprintf('%d',Kl);
fileC='_patches_';
fileD=sprintf('%d',n_pat);
fileE='_normal';
fileF='_resistant';
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
% plothours=plothours*tsph+1;

spaceplot=linspace(1,n_pat,n_pat);
% 

% figure(1);
% clf;
% subplot(2,2,1)
% for i=1:numplots
%     time=plothours(i);
%     ynorm=x_normal(time,:);
%     semilogy(spaceplot,ynorm,'Linewidth',2)
%     title('Sensitive Commensal Bacteria')
%     xlabel('Patches'); ylabel('Bacterial density')
% %     legend('t=0h', 't=6h', 't=12h', 't=24h', 't=48h', 't=72h', 't=96h', 't=720h')
%     legend('t=0h', 't=3d', 't=5d', 't=10d', 't=15d', 't=20d', 't=25d', 't=30d')
%     set(gca, 'Linewidth',2,'Fontsize',18)
%     hold on
% end
% 
% 
% subplot(2,2,2)
% for i=1:numplots
%     time=plothours(i);
%     yres=x_res(time,:);
%     semilogy(spaceplot,yres,'Linewidth',2)
%     title(' Resistant Commensal Bacteria')
%     xlabel('Patches');  ylabel('Bacterial density')
% %     legend('t=0h', 't=6h', 't=12h', 't=24h', 't=48h', 't=72h', 't=96h', 't=720h')
%     legend('t=0h', 't=3d', 't=5d', 't=10d', 't=15d', 't=20d', 't=25d', 't=30d')
%     set(gca, 'Linewidth',2,'Fontsize',18)
%     hold on
%     
% end
% 
% 
% subplot(2,2,3)
% for i=1:numplots
%     time=plothours(i);
%     ycommensal=x_res(time,:)+x_normal(time,:);
%     semilogy(spaceplot,ycommensal,'Linewidth',2)
%     title('All Commensal Bacteria')
%     xlabel('Patches'); ylabel('Bacterial density')
% %     legend('t=0h', 't=6h', 't=12h', 't=24h', 't=48h', 't=72h', 't=96h', 't=720h')
%     legend('t=0h', 't=3d', 't=5d', 't=10d', 't=15d', 't=20d', 't=25d', 't=30d')
%     set(gca, 'Linewidth',2,'Fontsize',18)
%     hold on
%     
% end
% 
% subplot(2,2,4)
% 
% for i=1:numplots
%     time=plothours(i);
%     yinfl=x_infl(time,:);
%     semilogy(spaceplot,yinfl,'Linewidth',2)
%     title('Inflammatory Bacteria')
%     xlabel('Patches'); ylabel('Bacterial density')
% %      legend('t=0h', 't=6h', 't=12h', 't=24h', 't=48h', 't=72h', 't=96h', 't=720h')
%     legend('t=0h', 't=3d', 't=5d', 't=10d', 't=15d', 't=20d', 't=25d', 't=30d')
%     set(gca, 'Linewidth',2,'Fontsize',18)
%     hold on
% end
% 
% %%%inflammatory response
% figure(2);
% clf;
% for i=1:numplots
%     time=plothours(i);
%     yil=Il(time,:);
%     plot(spaceplot,yil,'Linewidth',2)
%     title('Inflammatory response')
%     xlabel('Patches'); ylabel('Inflammation level')
% %     legend('t=0h', 't=6h', 't=12h', 't=24h', 't=48h', 't=72h', 't=96h', 't=720h')
%     legend('t=0h', 't=3d', 't=5d', 't=10d', 't=15d', 't=20d', 't=25d', 't=30d')
%     set(gca, 'Linewidth',2,'Fontsize',18)
%     hold on
%     
% end



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
title('Commensal Bacteria')
set(gca, 'Linewidth',2,'Fontsize',18)

subplot(2,2,2)
pcolor(spaceplot,t/240,log10(x_res));
shading interp;hold on;
[M,c]=contour(spaceplot,t/240,log10(x_res),10);
c.LineColor = 'k';
colormap parula
h=colorbar;
caxis([0 12]);
set(get(h,'label'),'string','Log_{10} Bacterial density');
xlabel('Patch'); ylabel('Time (days)');
title('Resistant Commensal Bacteria')
set(gca, 'Linewidth',2,'Fontsize',18)

subplot(2,2,3)
pcolor(spaceplot,t/240,log10(x_normal+x_res));
shading interp;hold on;
[M,c]=contour(spaceplot,t/240,log10(x_normal+x_res),10);
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
set(get(h,'label'),'string','Inflammation level');
xlabel('Patch'); ylabel('Time (days)');
title('Inflammatory response')
set(gca, 'Linewidth',2,'Fontsize',18)





%%%Violin plots in Fig 8

    
plothours=[0 72 120 240 1440 ];
plothours=plothours*tsph+1;

time=plothours(2);
ynorm=log10(x_normal(time,:));
ynorm(isinf(ynorm)) = 0; 
yres=log10(x_res(time,:));
yres(isinf(yres)) = 0; 
yinfl=log10(x_infl(time,:));
yinfl(isinf(yinfl)) = 0; 

time1=plothours(3);
ynorm1=log10(x_normal(time1,:));
ynorm1(isinf(ynorm1)) = 0; 
yres1=log10(x_res(time1,:));
yres1(isinf(yres1)) = 0; 
yinfl1=log10(x_infl(time1,:));
yinfl1(isinf(yinfl1)) = 0; 

time2=plothours(4);
ynorm2=log10(x_normal(time2,:));
ynorm2(isinf(ynorm2)) = 0; 
yres2=log10(x_res(time2,:));
yres2(isinf(yres2)) = 0; 
yinfl2=log10(x_infl(time2,:));
yinfl2(isinf(yinfl2)) = 0; 

time3=plothours(5);
ynorm3=log10(x_normal(time3,:));
ynorm3(isinf(ynorm3)) = 0; 
yres3=log10(x_res(time3,:));
yres3(isinf(yres3)) = 0; 
yinfl3=log10(x_infl(time3,:));
yinfl3(isinf(yinfl3)) = 0; 

X=[ynorm yres yinfl]';X1=[ynorm1 yres1 yinfl1]';X2=[ynorm2 yres2 yinfl2]';X3=[ynorm3 yres3 yinfl3]';

Y1=repmat('Sensitive',length(ynorm),1);size(Y1)
Y1 = cellstr(Y1);
Y2=repmat('Resistant',length(yres),1);size(Y2)
Y2 = cellstr(Y2);
Y3=repmat('Salmonella',length(yinfl),1);size(Y3)
Y3 = cellstr(Y3);
Bacteria=[Y1;Y2;Y3];


figure(11);
clf;
grouporder={'Sensitive', 'Resistant','Salmonella'};
subplot(2,2,1)
vs = violinplot(X, Bacteria,'GroupOrder',grouporder,'Width',0.2,'ViolinAlpha',0.4);
ylabel('Log_{10} Bacterial density');
title('Sensitive, resistant commensals and Inflammatory bacteria(Salmonella) at day 3')

subplot(2,2,2)
vs1 = violinplot(X1,Bacteria,'GroupOrder',grouporder,'Width',0.2,'ViolinAlpha',0.4);
ylabel('Log_{10} Bacterial density');
title('Sensitive, resistant commensals and Inflammatory bacteria(Salmonella) at day 5')

subplot(2,2,3)
vs2 = violinplot(X2,Bacteria,'GroupOrder',grouporder,'Width',0.2,'ViolinAlpha',0.4);
ylabel('Log_{10} Bacterial density');
 title('Sensitive, resistant commensals and Inflammatory bacteria(Salmonella) at day 10')

subplot(2,2,4)
vs3 = violinplot(X3, Bacteria,'GroupOrder',grouporder,'Width',0.2,'ViolinAlpha',0.4 );
ylabel('Log_{10} Bacterial density');
title('Sensitive, resistant commensals and Inflammatory bacteria(Salmonella) at day 60')
end

