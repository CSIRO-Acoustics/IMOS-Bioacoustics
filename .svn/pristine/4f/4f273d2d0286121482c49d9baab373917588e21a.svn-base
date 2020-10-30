

[File, Path] = uigetfile(pwd,'Navigate to the required data_file_list.txt file','*.txt');
dataFile = fullfile(Path, File)
DestinationDir = uigetdir('Q:\Processed_data', 'Navigate to the destination directory');

%% copy files from the data_file_list.txt file to a nominated destination folder

fid = fopen(dataFile,'r');

while ~feof(fid)
      sourceFile = fgetl(fid);
      %keyboard
      srtInd = max(strfind(sourceFile,'\')+1);
      fileName = sourceFile(srtInd:end);
      destinfile = [DestinationDir '\' fileName];
      copyfile(sourceFile, destinfile);
      fprintf('Copied %s to %s\n',sourceFile,destinfile)
end

close all
clear all


%%  Old code from copy4sure.mat which has some nice parts: looping till 
%   identified file has copied and reports how many files have failed to copy and how 
%   many have been successfully copied.
    

    % SourceDir = uigetdir(pwd, 'Navigate to source files');
    % cd(SourceDir);


    % if ~isdir(DestinationDir)
    %     mkdir(DestinationDir)
    % end 
    % 
    % [m,n]=size(filelist);
    % SUCCESS = 0; num_failures = 0;
    % 
    % 
    % for i=3:m    
    %     %SourceFile = [SourceDir '\' strtrim(Source_files(i,:))];
    %     SourceFile = filelist(i);
    %     DestinationFile = [DestinationDir '\' filelist(i)];
    %     if ~isdir(SourceFile)        
    %         while isequal(SUCCESS,0)
    %             [SUCCESS,MESSAGE,MESSAGEID] = copyfile(SourceFile,DestinationFile);           
    %             if isequal(SUCCESS,1)
    %                 fprintf('Copied %s to %s\n',SourceFile,DestinationFile)
    %             else
    %                 fprintf('Failed to copy %s to %s, try again\n',SourceFile,DestinationFile)
    %                 num_failures = num_failures +1;
    %             end
    %         end
    %     end
    %     SUCCESS = 0;
    % end