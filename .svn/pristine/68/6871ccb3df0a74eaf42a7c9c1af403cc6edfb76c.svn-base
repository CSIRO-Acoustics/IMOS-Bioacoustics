function make_datafile_list(datadir)
% 
% small function to generate a textfile that contains the full path + name
% of all raw files in the data directory
%
% Tim ryan 10/05/2017

cd(datadir);
rawfiles = ls('*.raw');
fid = fopen('Raw_data_files.txt','w+');
for i=1:size(rawfiles,1)
    fprintf(fid,'%s\n', [pwd '\' rawfiles(i,:)]);
end
fclose('all');
