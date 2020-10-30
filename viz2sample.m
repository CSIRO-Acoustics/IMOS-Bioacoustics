function sample_data = viz2sample(data)
% VIZ2SAMPLE converts a data structure returned by viz_sv into IMOS toolbox
% sample_data format.

D=1;
T=2;
C=3;

sample_data=struct();

ld = length(data.depth);
lt = length(data.time);
lc = 1;

sample_data.dimensions{D}.name = 'DEPTH';
sample_data.dimensions{D}.data = data.depth;
sample_data.dimensions{T}.name = 'TIME';
sample_data.dimensions{T}.data = data.time;
sample_data.dimensions{C}.name = 'Channel';
sample_data.dimensions{C}.data = {};

if isfield(data, 'channels')
    lc = length(data.channels);
    sample_data.dimensions{C}.data = data.channels;
    data=rmfield(data,'channels');
end

sample_data.variables{1}.name = 'LONGITUDE';
sample_data.variables{1}.dimensions = T;
sample_data.variables{1}.data = data.longitude;

sample_data.variables{2}.name = 'LATITUDE';
sample_data.variables{2}.dimensions = T;
sample_data.variables{2}.data = data.latitude;

data=rmfield(data,{'depth','time','latitude','longitude'});

fields = fieldnames(data);
for fld = fields'
    field = fld{1};
    if isscalar(data.(field)) || ischar(data.(field))
        sample_data.(field) = data.(field);
    else
        sample_data.variables{end+1}.name = field;
        sample_data.variables{end}.data = data.(field);
        sz = size(data.(field));
        sz(sz == 1) = [];
        for i = 1:length(sz)
            if sz(i) == lt
                sample_data.variables{end}.dimensions(i) = T;
            elseif sz(i) == ld
                sample_data.variables{end}.dimensions(i) = D;
            elseif sz(i) == lc
                sample_data.variables{end}.dimensions(i) = C;
            else
                sample_data.dimensions{end+1}.name = [field num2str(i)];
                sample_data.dimensions{end}.data = 1:sz(i);
                sample_data.variables{end}.dimensions(i) = length(sample_data.dimensions);
            end
        end
    end
end
                
                