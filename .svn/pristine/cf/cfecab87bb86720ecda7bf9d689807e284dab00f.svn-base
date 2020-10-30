function stripped = cut_depth( data, depth )
% cut_depth trims a viz_sv data structure to the depth specified.

stripped = data;

include = data.depth <= depth;
fields = fieldnames(data);
for f=1:length(fields)
    if size(data.(fields{f}), 1) == length(include)
        stripped.(fields{f}) = data.(fields{f})(include,:,:);
    end
end



