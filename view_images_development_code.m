%%
% This is some rough and ready code to look at possibilities of viewing
% echoview jpg images 
% Further development is planned. 
% Tim Ryan 18/10/2019
close all
%%
image_folder = 'Q:\Pending_registration\Year2019\Tokatu\Tokatu_received_20191001\images';
filenames = dir(fullfile(image_folder, '*.jpg'));  % read all images with specified extention, its jpg in our case
total_images = numel(filenames);    % count total number of photos present in that folder
%for n = 1:total_images
f1 = figure(1) ;
f2 = figure(2) ;
f1.Position=[525 95 1455 893]
%get(0,'MonitorPositions')

n=1; 
%for n = 1:100   
%manual = 'False';
manual = 'True';
set(gcf,'currentchar',' ')         % set a dummy character

while n<100     
    full_name= fullfile(image_folder, filenames(n).name);         % it will specify images names with full path and extension
    f=full_name;
    raw_image = imread(f);                 % Read images            
    cleaned_image = imread(f);             % Read images            
    figure(1);
    %subplot(211)
    imshow(raw_image);                     % Show all images   
    %subplot(212)
    figure(2);    
    imshow(cleaned_image)
    f1.Position=[-1525 95 1455 893];
    f2.Position=[525 95 1455 893];
    %pause
    warning off
    WinOnTop(f1,true);
    %WinOnTop(f2,true);
    warning on
    if isequal(manual,'True')
        k = waitforbuttonpress;
        % 28 leftarrow
        % 29 rightarrow
        % 30 uparrow
        % 31 downarrow
        value = double(get(gcf,'CurrentCharacter'));
        if isequal(value,28)
            n=n-1;
        else
            n=n+1;
        end
        fprintf('value is %5d, n is %d\n',value,n);        
        clf
    else
        pause(0.01);
        n=n+1;    
        halt = get(gcf,'currentchar')==' ' ;             
        click_type=get(f1,'SelectionType');
        if isequal(click_type,'alt')
            fprintf('type dbcont to continue\n');
            keyboard            
            click_type = 'normal';
        end
    end
end
 
%%
while get(gcf,'currentchar')==' '  % which gets changed when key is pressed
   do_stuff()
end