%%%%Code used to generate Fig 7, 8 & 9
function s_initialconditions_resitant
% Read in no infection results, and produce an initial conditions file for
% with infection model.

% Set parameters
tsph=10;        % time steps per hour
tth=1440;         % total time in hours
tt=tth*tsph;
Ks = 1e9;     % Carrying capacity in small intestine
Kl = 1e12;    % Carrying capacity in large intestine
n_pat=60;
results_st=240;

% Set up filename
fileA='Bacteria_n_only_Ks_';
fileB=sprintf('%d',Ks);
fileBX='_Kl_';
fileBY=sprintf('%d',Kl);
fileC='_patches_';
fileD=sprintf('%d',n_pat);
fileE='_normal';
fileF='_resistant'; %%added
fileX='_end';
fileY=sprintf('%d',tt);
fileG='.mat';

fileM1={fileA,fileB,fileBX,fileBY,fileC,fileD,fileX,fileY,fileE,fileG};
filename1=strjoin(fileM1,'');

%%%added
fileM2={fileA,fileB,fileBX,fileBY,fileC,fileD,fileX,fileY,fileF,fileG};
filename2=strjoin(fileM2,'');

load(filename1)
load(filename2)

InCon = zeros(1,n_pat);
InConR = zeros(1,n_pat);

for time=results_st:tt
    for i=1:n_pat
        InCon(i)=InCon(i)+x_normal(time+1,i);
        InConR(i)=InConR(i)+x_res(time+1,i);
    end
end

InCon=round(InCon./(tt-results_st+1));
InConR=round(InConR./(tt-results_st+1));

file1A='InitialConditions_Ks_';
fileG='.mat';
file1M={file1A,fileB,fileBX,fileBY,fileC,fileD,fileG};
filename1a=strjoin(file1M,'');
save(filename1a, 'InCon');

file1AR='InitialConditions_R_Ks_';
fileG='.mat';
file1MR={file1AR,fileB,fileBX,fileBY,fileC,fileD,fileG};
filename1aR=strjoin(file1MR,'');
save(filename1aR, 'InConR');

end