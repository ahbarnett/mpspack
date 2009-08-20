function s = merge(s1,s2)
% MERGE - Merges structures together
%
%  s = MERGE(s1,s2) returns a struct that contains all fields 
%  from the structs s1 and s2. If s1 and s2 have a field with 
%  the same name s1 has preference.

s=s1;
names=fieldnames(s2);
for j=1:length(names),
    field=names{j};
    if ~isfield(s,field),
        s=setfield(s,field,getfield(s2,field));
    end
end
