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
filenames=ls('*.raw');
output_file = 'test.txt';
%filenames = 'SXO-D20090715-T132350.raw'

dir_write=dir_read;
%dir_write='Z:\temp\matlab\read_write_ek60\' % directory to write from
fid3 = fopen([dir_write output_file],'w');
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
%fname=strcat(dir_read,file_name);
%fname2=strcat(dir_write,file_name(1:21),num2str(split_on),'.raw'); % file to write - 1st 24 letters 
[m,n]=size(filenames);
for i=1:m   
    filename = strtrim(filenames(i,:));
    fid = fopen(filename,'r');
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
        try
            length_header(2)=length; % length of header
        catch
            keyboard
        end

        % Read NMEA, Annotation, or Sample datagram
        while (1)
            length = fread(fid,1,'int32');
            if (feof(fid)) % if end of file
                break
            end
            dgheader = readdgheader(fid);
            dgheader.datagramtype;
            switch (dgheader.datagramtype)
            case 'NME0' % NMEA datagram
                text = readtextdata(fid,length-headerlength);                  
                nmea_time = datestr((dgheader.datetime/10^7)/(24*60*60) + datenum(1601,1,1));
                true_time = ((dgheader.datetime/10^7)/(24*60*60) + datenum(1601,1,1));
                millisec = round(1000*(true_time - datenum(nmea_time))*24*60*60);
    %           if text.text(1:6) == '$PAMOT' % we have a pr telegram
               if text.text(1:6) == '$PADPT' % we have a depth telegram        
                    [A, Platform_depth, C] = strread(text.text,'%s%f%s','delimiter',',');
                    nmea_time  = [datestr(datenum(nmea_time),'yyyymmdd HHMMSS') sprintf('%03d', (millisec)) '0' ];               
                    try
                        fprintf(fid3, '%s %f 3\n', nmea_time,Platform_depth);
                    catch
                        keyboard
                    end
               end 

             case 'TAG0' % Annotation datagram
                text = readtextdata(fid,length-headerlength)
                writetextdata(fid2,text); % ADR write
             case 'RAW0' % Sample datagram
                sampledata = readsampledataraw(fid);
                 %writesampledata(fid2,sampledata);
                pingno = pingno +1;
             otherwise
                fprintf('unknown datagram in file\n');
                fprintf('skip to the next file\n');
                break
%                    error(strcat('Unknown datagram ''',dgheader.datagramtype,''' in file'));
            end
            length = fread(fid,1,'int32');
        end
        fclose(fid);
        fprintf('finished processing file %s\n',filename)
    end
end
disp('Finished reading file - Happy Pinging and Singing');
    
%% now insert the header into the file
fid = fopen(output_file,'r')
j=0;p=1;;
while p > -1
    p=fgetl(fid)
    j=j+1;
end
fclose(fid)
fid2 = fopen([output_file(1:end-4) '.evl'],'w');
fprintf(fid2, 'EVBD 3 4.50.46.12102\n');
fprintf(fid2, '%s\n',num2str(j-1));
fclose(fid2);

% now cat these together
a = [output_file(1:end-4) '.evl']
%a = output_file;
b = 'test.txt'
system(['for %f in ("' b '") do type "%f" >> "' a '"']);

%%



    
