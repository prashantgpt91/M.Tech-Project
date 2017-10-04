% make directories
% make directory structures, for concept detection and duplicate detection

RootDir = GetRootPath;

ret = mkdir(RootDir,'Corners');
ret = mkdir(RootDir,'Features');
ret = mkdir(RootDir,'Gargs');
ret = mkdir(RootDir,'Images');
ret = mkdir(RootDir,'Lists');
ret = mkdir(RootDir,'Results');
ret = mkdir([RootDir 'Results\'],'X');

ret = mkdir(RootDir,'Training');
ret = mkdir([RootDir 'Training\'],'X');
ret = mkdir([RootDir 'Training\'],'Parameters');
