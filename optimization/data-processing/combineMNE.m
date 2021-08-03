% script to combine tagged and all cells from Economo 2018 Nature paper

data_pth = '/Users/Munib/Documents/Economo-Lab/subspace-id/data/elsayed/mike';
% load all cell data
load(fullfile(data_pth,'alldat.mat'));
% load pt upper and lower cell data
load(fullfile(data_pth,'ptlow.mat'));
load(fullfile(data_pth,'ptup.mat'));

tag.psth = [tag.psth;ptlow;ptup];
tag.ptlow_ix = 1208:1276;
tag.ptup_ix = 1277:1337;

save(fullfile(data_pth,'alldat_tagged.mat'),'tag','-v7.3');

