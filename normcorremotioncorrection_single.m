function normcorremotioncorrection_single(name,tempname)
%Modified for a use in single channel tifstack
addpath(genpath(fullfile('/gpfs/scratch/karadm01','NoRMCorre-master')));

[pathname,filename,ext] = fileparts(name);
tempPath = fullfile(pathname,tempname);

tic; Y = read_file(name); toc;
G = double(Y);      % convert to double precision 
%[G,R] = deinterleave(Y); %deinterleave into green and red channels
T = size(G,ndims(G));
template = read_file(tempPath);
% [template,~] = deinterleave(template);

options_nonrigid = NoRMCorreSetParms('d1',size(G,1),'d2',size(G,2),...
    'grid_size',[32,32],'mot_uf',4,'bin_width',10,'max_shift',15,...
    'max_dev',10,'us_fac',10,'upd_template',false);
tic; 
[M2,shifts2,template2] = normcorre_batch(G,options_nonrigid,template); 
toc
M2 = uint16(M2);

if ~exist(fullfile(pathname,'aligned'), 'dir')
   mkdir(fullfile(pathname,'aligned'))
end

saveastiff(M2,fullfile(pathname,'aligned',[filename ext]));

function [A, B] = deinterleave(I)
    A = I(:,:,1:2:end);
    B = I(:,:,2:2:end);
end
end