function sots_reprocess(ncfilelist)

if nargin < 1 || isempty(ncfilelist)
    [file,path] = uigetfile('*.nc', 'multiselect', 'on');
    ncfilelist = fullfile(path,file);
end

if ischar(ncfilelist) && exist(ncfilelist,'file') == 2 && ~strcmp('.nc', ncfilelist(end-2:end))
    fid = fopen(ncfilelist);
    ncfilelist = textscan(fid,'%s');
    ncfilelist = ncfilelist{1};
end
if ischar(ncfilelist)
    ncfilelist = {ncfilelist};
end

settings = load('sots_reprocess.mat');
for ncfile = ncfilelist'
    reprocessev(settings,ncfile{1});
end
