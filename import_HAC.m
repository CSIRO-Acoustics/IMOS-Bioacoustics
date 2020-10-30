function import_HAC(EvApp, directory, ev_filelist, FilesetName, VarName)
% Import the requested variable from the specified worksheets to HAC files.
%
% Inputs:
%   EvApp       a handle to a COM object of an EchoviewApplication.
%   directory   name of directory containing HAC files 'HAC_outputs'
%   ev_filelist list of worksheets to import into    
%   FilesetName EchoView Fileset name 'HAC1' - 'HAC5'
%   VarName     name of variable to import 'HAC1_exp' - 'HAC5_exp'
%

    % -- Setup HAC output directory
    [path, ~, ~] = fileparts(ev_filelist{1});
    [path, ~, ~] = fileparts(path);
    
    HAC_Output_dir = fullfile(path, directory);
    
    for i=1:length(ev_filelist)
        [~, name, ext] = fileparts(ev_filelist{i});
        hac_file_to_add = fullfile(HAC_Output_dir, [ name ' ' VarName '_000000.hac']);
        fprintf('\nImporting HAC file %d of %d:\n%s \ninto\nev fileset:\n%s\n\n', ...
            i, length(ev_filelist), hac_file_to_add, [name ext]);
        
        EvFile = EvApp.OpenFile(ev_filelist{i});        
        
        try
            if ~EvFile.Filesets.FindByName(FilesetName).DataFiles.Add(hac_file_to_add);
                fprintf('File not added - check the name and path is correct\n')
            end
        catch exception
            fprintf('Stopped at %s line %d because of:\n%s\n', ...
                exception.stack(1).file, exception.stack(1).line, exception.message);
            fprintf('type "dbcont" and press Enter to continue\n');
            keyboard           
        end
        EvFile.Save;
        EvFile.Close;
    end
end
