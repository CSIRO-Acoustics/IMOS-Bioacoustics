%read_write_raw.m
% Program to read in EK60 and ES60 Raw files and split them into smaller
% files
% to use, alter user inputs below and execute - this results in shorter files
% last part of fienames gives the order of files
%
% Alex De Robertis AFSC 2/18/03 - based on code by Lars Nonboe Andersen
% of simrad

clear all
close all

% USER INPUTS %%%%%%%%%%%%%%%%%%%  DO NOT EDIT BELOW THIS LINE %%%%%%%%
%dir_read='Z:\temp\matlab\read_write_ek60\' % directory to read from
dir_read = 'H:\Year2009_Data\SaxonOnward_2009\SXO200901\OP2\Acoustics\'
file_name='SXO-D20090715-T130857.raw' %file to read
dir_write=dir_read
%dir_write='Z:\temp\matlab\read_write_ek60\' % directory to write from
fid3 = fopen([dir_write 'test.txt'],'w');

%dir_read='C:\Data\ES60_gulf\' % directory to read from
%file_name='L0001-D20030716-T104119-ES60.raw' %file to read
%dir_write='C:\Data\ES60_gulf\' % directory to write from
%num_pings=360000;  % incremented number of pings to split files
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
 
% writeheader;  % write the header material to new file
    fprintf(fid3, 'date_time, pitch,roll\n')
    % Read NMEA, Annotation, or Sample datagram
    while (1)
        
       % if ((pingno-pinglast)==num_pings);  % check if should split file
       %     pinglast=pingno;
       %     split_on=split_on+1;
       %     fname2=strcat(dir_write,file_name(1:21),num2str(split_on),'.raw'); % file to write - 1st 24 letters 
       %     fclose(fid2);
       %     fid2=fopen(fname2,'w'); % open up file to write to 
       %     disp(fname2)
       %     writeheader;  % write the header material to new file
       % end
        
        length = fread(fid,1,'int32');
        fwrite(fid2,length,'int32'); % ADR write
        
        if (feof(fid)) % if end of file
            break
        end
        
        dgheader = readdgheader(fid);
%        writedgheader(fid2,dgheader); % ADR write
        dgheader.datagramtype;
        switch (dgheader.datagramtype)
        case 'NME0' % NMEA datagram
            text = readtextdata(fid,length-headerlength);                  
            nmea_time = datestr((dgheader.datetime/10^7)/(24*60*60) + datenum(1601,1,1));
            true_time = ((dgheader.datetime/10^7)/(24*60*60) + datenum(1601,1,1));
            millisec = round(1000*(true_time - datenum(nmea_time))*24*60*60);
           %  writetextdata(fid2,text); % ADR write
%           if text.text(1:6) == '$PAMOT' % we have a pr telegram
           if text.text(1:6) == '$PADPT' % we have a pr telegram        
               keyboard
               [A B pitch D roll]=strread(text.text,'%s%s%f%s%f','delimiter',',');
                fprintf(fid3, '%s.%0.3d, %.2f, %.2f\n', nmea_time, millisec, pitch, roll);
           end 
           
         case 'TAG0' % Annotation datagram
             text = readtextdata(fid,length-headerlength)
            writetextdata(fid2,text); % ADR write
 
         case 'RAW0' % Sample datagram
            sampledata = readsampledataraw(fid);
             %writesampledata(fid2,sampledata);
             pingno = pingno +1;
        otherwise
            error(strcat('Unknown datagram ''',dgheader.datagramtype,''' in file'));
        end
        length = fread(fid,1,'int32');
%         fwrite(fid2,length,'int32'); % ADR write

    end
    fclose(fid);
    fclose(fid2);
end

disp('Finished reading file - Happy Pinging and Singing');
    
    
