% test qp bounded strip layer-potential basis. Barnett 3/12/10
clear all classes; v=1; % verbosity
thi = -pi/5; om = 10; %7.76644415490187; %(t.a=1) % inc ang, overall freq
M = 35; Mt = 25;
d = 1; yB = -1; yT = 1;       % strip & box geom
B = segment(Mt, 1i*yB+[-d/2,d/2], 'g'); T = segment(Mt, 1i*yT+[d/2,-d/2], 'g');
t = boundedqpstrip([B T], om, struct('M', M));
kvec = om*exp(1i*thi); ui = @(x) exp(1i*real(conj(kvec) .* x)); 
a = exp(1i*real(conj(kvec) * d)); t.setbloch(a);
t.addqpbstlayerpots;
co = kron([0;1], ones(t.bas{1}.Nf,1));

if v, dx = 0.025; ep=1e-3; x = -1.5+ep:dx:1.5; y = -2:dx:2;% plot grid (shifted)
[xx yy] = meshgrid(x,y); p = pointset(xx(:)+1i*yy(:));
ug = reshape(t.evalbases(p) * co, size(xx));
figure; imagesc(x, y, real(ug)); set(gca,'ydir','normal'); axis tight equal;
caxis([-1 1]); colormap(jet(256)); colorbar; t.showbasesgeom;
end

Q = t.evalbasesdiscrep;
%figure; imagesc(real(Q)); axis equal tight; colormap(jet(256)); colorbar;

