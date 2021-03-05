function s_initialconditions
% Read in no infection results, and produce an initial conditions file for
% with infection model.

% Set parameters
tsph=10;        % time steps per hour
tth=720;         % total time in hours
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
fileX='_end';
fileY=sprintf('%d',tt);
fileG='.mat';

fileM1={fileA,fileB,fileBX,fileBY,fileC,fileD,fileX,fileY,fileE,fileG};
filename1=strjoin(fileM1,'');

load(filename1)

InCon = zeros(1,n_pat);

for time=results_st:tt
    for i=1:n_pat
        InCon(i)=InCon(i)+x_normal(time+1,i);
    end
end

InCon=round(InCon./(tt-results_st+1));

file1A='InitialConditions_Ks_';
fileG='.mat';
file1M={file1A,fileB,fileBX,fileBY,fileC,fileD,fileG};
filename1a=strjoin(file1M,'');
save(filename1a, 'InCon');

end