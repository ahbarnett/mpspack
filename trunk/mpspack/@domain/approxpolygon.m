% APPROXPOLYGON - extracts col vec of approx vertices from closed seg array
%
function v = approxpolygon(s, pm)
v = [];
for j=1:length(s)
  if pm(j)==1
    v = [v; s(j).approxv(1:end-1)];
  else
     v = [v; s(j).approxv(end:-1:2)];
  end
end
