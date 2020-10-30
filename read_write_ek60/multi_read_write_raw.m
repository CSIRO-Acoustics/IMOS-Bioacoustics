%mliti read_write_raw.m
% Program to read in EK60 and ES60 Raw files and split them into smaller
% files
% to use, alter user inputs below and execute - this results in shorter files
% last part of fienames gives the order of files
%
% Alex De Robertis AFSC 2/18/03 - based on code by Lars Nonboe Andersen
% of simrad

clear all
close all



%files=dir ('C:\Data\ES60_gulf\raw_data\*.raw')  % directory command to look for files to process

files=dir ('Z:\temp\matlab\read_write_ek60\*.raw')  % directory command to look for files to process


for i=1:length(files)  % loop for all files in directory
    % USER INPUTS %%%%%%%%%%%%%%%%%%%  
    dir_read='Z:\temp\matlab\read_write_ek60\' % directory to read from
    file_name=files(i).name %file to read
    dir_write='Z:\temp\matlab\read_write_ek60\' % directory to write new data to
    num_pings=7200;  % incremented number of pings to split files
    %%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    % constants
    headerlength = 12; % Bytes in datagram header
    pingno = 1;% ping currently on
    pinglast=1;% counter for ping last used to split file
    split_on=1; % file number
    
    % open files
    fname=strcat(dir_read,file_name);

    fname2=strcat(dir_write,file_name(1:21),num2str(split_on),'.raw'); % file to write - 1st 24 letters 
    fid = fopen(fname,'r');
    fid2=fopen(fname2,'w'); % open up file to write to write to
    
    if (fid==-1)
        error('Could not open file');
    else
        % Read configuration datagram
        length = fread(fid,1,'int32');
        length_header(1)=length; % length of header
        dgheader = readdgheader(fid);
        configheader = readconfigheader(fid);
        for i=1:configheader.transducercount;
            configtransducer(i) = readconfigtransducer(fid);
        end
        config = struct('header',configheader,'transducer',configtransducer);
        length = fread(fid,1,'int32');
        length_header(2)=length; % length of header
        
        keyboard
        writeheader;  % write the header material to new file
        
        % Read NMEA, Annotation, or Sample datagram
        while (1)
            
            if ((pingno-pinglast)==num_pings);  % check if should split file
                pinglast=pingno;
                split_on=split_on+1;
                fname2=strcat(dir_write,file_name(1:24),num2str(split_on),'.raw'); % file to write - 1st 24 letters 
                fclose(fid2);
                fid2=fopen(fname2,'w'); % open up file to write to 
                disp(fname2)
                writeheader;  % write the header material to new file
            end
            
            length = fread(fid,1,'int32');
            fwrite(fid2,length,'int32'); % ADR write
            
            if (feof(fid)) % if end of file
                break
            end
            
            dgheader = readdgheader(fid);
            writedgheader(fid2,dgheader); % ADR write
            switch (dgheader.datagramtype)
                case 'NME0' % NMEA datagram
                    text = readtextdata(fid,length-headerlength);
                    writetextdata(fid2,text); % ADR write                 
                case 'TAG0' % Annotation datagram
                    text = readtextdata(fid,length-headerlength);
                    writetextdata(fid2,text); % ADR write
                case 'RAW0' % Sample datagram
                    sampledata = readsampledataraw(fid);
                    writesampledata(fid2,sampledata);% adr write
                    pingno = pingno +1;                   
                otherwise
                    error(strcat('Unknown datagram ''',dgheader.datagramtype,''' in file'));
            end
            length = fread(fid,1,'int32');
            fwrite(fid2,length,'int32'); % ADR write
            
        end
        fclose(fid);
        fclose(fid2);
    end
end % for major loop
disp('Finished reading file - Happy Pinging and Singing');


