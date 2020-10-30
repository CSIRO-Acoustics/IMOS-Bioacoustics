%writeheader.m
% write header material for ek60 data file


    fwrite(fid2,length_header(1),'int32'); % ADR write length
    writedgheader(fid2,dgheader); % ADR write dgheader
    writeconfigheader(fid2,configheader); % ADR write configheader
    for i=1:configheader.transducercount
        writeconfigtransducerraw(fid2,configtransducer);  % ADR write congfigtransducer
    end
    fwrite(fid2,length_header(2),'int32'); % ADR write