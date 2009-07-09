% MONOCHROME - make array of graphics handles go one color
%
function monochrome(h, c)
  hl = [findall(h,'Type', 'line'); findall(h,'Type', 'text'); ...
        findall(h, 'Type', 'marker')]; set(hl, 'color', c);
  hp = findall(h,'Type', 'patch'); set(hp, 'FaceColor', c, 'Edgecolor', c);
