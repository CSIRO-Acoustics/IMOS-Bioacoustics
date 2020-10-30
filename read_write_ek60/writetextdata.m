% Reading EK60 raw data file text data
% Simrad, Lars Nonboe Andersen, 21/12-01

function writetextdata(fid2,text);
% length in bytes

%textdata.text = char(fread(fid,length,'char')');
fwrite(fid2,text.text,'char');