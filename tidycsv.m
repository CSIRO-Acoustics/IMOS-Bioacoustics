    function tidycsv(datadir)
% utility to remove 'no data' rows of data in the csv files. This is a
% workaround for the fact that the Echoview investigator template is outputting no data
%
%
% Tim Ryan
% 10-Jan-2016
%
% 
csvs = ls([datadir '\*.csv']);
for i=1:size(csvs,1)
    fprintf('processing %s',csvs(i,:));    
    clean_csv(fullfile(datadir, strtrim(csvs(i,:))));
    fprintf('Done\n');
end

cleanup = questdlg('Would you like to clean up the old files?','Clean up old files','Yes','No','Yes');
if isequal(cleanup,'Yes')
    delete([datadir '\*.csv']);
    txtfiles = ls([datadir '\*.txt']);
    for i=1:size(txtfiles,1)
        oldfile = [datadir '\' txtfiles(i,:)];
        [ndir, nfile, newext] = fileparts([datadir '\' txtfiles(i,:)]);
        newfile = fullfile(ndir, [nfile, '.csv']);
        %keyboard
        movefile(oldfile,newfile);
    end
end




function clean_csv(csvfile)
    fid = fopen(csvfile,'r');
    [dir,file,ext] = fileparts(csvfile);
    file1 = [file '.txt'];
    csvfile1 = fullfile(dir,file1);
    fid1 = fopen(csvfile1,'w+');
    ln = fgetl(fid);
    ln = strrep(ln,ln(regexpi(ln, '[^\s\x{20}-\x{7e}]')),'');
    fprintf(fid1, '%s\n',ln);
    while ~feof(fid)
        ln = fgetl(fid);
        ln = strrep(ln,ln(regexpi(ln, '[^\s\x{20}-\x{7e}]')),'');
        k = strsplit(ln,',');
        Sv = str2num(k{4});           
        if ~(isequal(Sv,9999) | isequal(Sv,-999))
            fprintf(fid1,'%s\n',ln);
        end
    end
    fclose('all');



        