% Timo's attempt to handle the cavity...

M=190;
k=25;
N=190;
rmfs=7;


% Now define all 81 segments...

l1=2+6i-(1+2i); l2=4+6i-(5+2i);


s=segment.empty(89,0);

s(1)=segment(M,[1+2i 2+2i]);
s(2)=segment(M,[2+2i 4+2i]);
s(3)=segment(M,[4+2i 5+2i]);
s(4)=segment(M,1+2i+[0 l1/4]);
s(5)=segment(M,[2+2i 2+3i]);
s(6)=segment(M,[4+2i 4+3i]);
s(7)=segment(M,[5+2i 5+2i+l2/4]);
s(8)=segment(M,5+2i+[l2/4 l2/4*2]);
s(9)=segment(M,[4+3i 5+2i+l2/4]);
s(10)=segment(M,[3+3i 4+3i]);
s(11)=segment(M,[3+3i 3+4i]);
s(12)=segment(M,[2+3i 3+3i]);
s(13)=segment(M,[1+2i+l1/4 2+3i]);
s(14)=segment(M,1+2i+[l1/4 l1/4*2]);
s(15)=segment(M,1+2i+[l1/4*2 l1/4*3]);
s(16)=segment(M,[1+2i+l1/4*2 3+4i]);
s(17)=segment(M,[3+4i 3+5i]);
s(18)=segment(M,[3+4i 5+2i+l2/4*2]);
s(19)=segment(M,[5+2i+l2/4*2 5+2i+l2/4*3]);
s(20)=segment(M,5+2i+[l2/4*3 l2/4*3.5]);
s(21)=segment(M,[3+5i 5+2i+l2/4*3]);
s(22)=segment(M,[3+5i 3+5.5i]);
s(23)=segment(M,[1+2i+l1/4*3 3+5i]);
s(24)=segment(M,1+2i+[l1/4*3 l1/4*3.5]);
s(25)=segment(M,1+2i+[l1/4*3.5 l1]);
s(26)=segment(M,[1+2i+l1/4*3.5 3+5.5i]);
s(27)=segment(M,[3+5.5i 5+2i+l2/4*3.5]);
s(28)=segment(M,5+2i+[l2/4*3.5 l2]);
s(29)=segment(M,[1.5+6i 2+6i]);
s(30)=segment(M,[3+5.5i 3+6.5i]);
s(31)=segment(M,[4+6i 4.5+6i]);
s(32)=segment(M,[4.5+6i 4.5+6.5i]);
s(33)=segment(M,[0.5+6i 0.5+6.5i]);
s(34)=segment(M,[0.5+6.5i 1.5+6.5i]);
s(35)=segment(M,[1.5+6.5i 3+6.5i]);
s(36)=segment(M,[3+6.5i 4.5+6.5i]);
s(37)=segment(M,[4.5+6.5i 5.5+6.5i]);
s(38)=segment(M,[5.5+6i 5.5+6.5i]);
s(39)=segment(M,[0.5+6i 1.5+6i]);
s(40)=segment(M,[4.5+6i 5.5+6i]);
s(41)=segment(M,[1.5+6i 1.5+6.5i]);
s(42)=segment(M,[0.5+6.5i 0.5+7.5i]);
s(43)=segment(M,[1.5+6.5i 1.5+7.5i]);
s(44)=segment(M,[3+6.5i 3+7.5i]);
s(45)=segment(M,[4.5+6.5i 4.5+7.5i]);
s(46)=segment(M,[5.5+6.5i 5.5+7.5i]);
s(47)=segment(M,[0.5+7.5i 1.5+7.5i]);
s(48)=segment(M,[1.5+7.5i 3+7.5i]);
s(49)=segment(M,[3+7.5i 4.5+7.5i]);
s(50)=segment(M,[4.5+7.5i 5.5+7.5i]);
s(51)=segment(M,[0+6i 0.5+6i]);
s(52)=segment(M,[0+5.5i 0+6i]);
s(53)=segment(M,[-1.5+5.5i 0+5.5i]);
s(54)=segment(M,[-1.5+5.5i -1.5+7.5i]);
s(55)=segment(M,[-1.5+7.5i 0.5+7.5i]);
s(56)=segment(M,[6+5.5i 6+6i]);
s(57)=segment(M,[6+5.5i 7.5+5.5i]);
s(58)=segment(M,[7.5+5.5i 7.5+7.5i]);
s(59)=segment(M,[5.5+7.5i 7.5+7.5i]);
s(60)=segment(M,[-1.5+4.5i -1.5+5.5i]);
s(61)=segment(M,[0+4.5i 0+5.5i]);
s(62)=segment(M,[-1.5+4.5i 0+4.5i]);
s(63)=segment(M,[-1.5+3i -1.5+4.5i]);
s(64)=segment(M,[-1.5+3i 0+3i]);
s(65)=segment(M,[0+3i 0+4.5i]);
s(66)=segment(3*M,[-1.5-1.5i -1.5+3i]);
s(67)=segment(2*M,[0 3i]);
s(68)=segment(3*M,[-1.5-1.5i 3-1.5i]);
s(69)=segment(2*M,[0 3]);
s(70)=segment(3*M,[3-1.5i 7.5-1.5i]);
s(71)=segment(3*M,[7.5-1.5i 7.5+3i]);
s(72)=segment(2*M,[3 6]);
s(73)=segment(2*M,[6 6+3i]);
s(74)=segment(M,[6+3i 6+4.5i]);
s(75)=segment(M,[6+3i 7.5+3i]);
s(76)=segment(M,[7.5+3i 7.5+4.5i]);
s(77)=segment(M,[6+4.5i 6+5.5i]);
s(78)=segment(M,[6+4.5i 7.5+4.5i]);
s(79)=segment(M,[7.5+4.5i  7.5+5.5i]);
s(80)=segment(M,[3-1.5i 3]);
s(81)=segment(M,[5.5+6i 6+6i]);
s(82)=segment(M,3+3i+[4.5*sqrt(2)*exp(1i*pi/4) rmfs*exp(1i*pi/4)]);
s(83)=segment(M,3+3i+[4.5*sqrt(2)*exp(3i*pi/4) rmfs*exp(3i*pi/4)]);
s(84)=segment(M,3+3i+[4.5*sqrt(2)*exp(5i*pi/4) rmfs*exp(5i*pi/4)]);
s(85)=segment(M,3+3i+[4.5*sqrt(2)*exp(-1i*pi/4) rmfs*exp(-1i*pi/4)]);
s(86)=segment(3*M,[3+3i rmfs pi/4 3*pi/4]);
s(87)=segment(3*M,[3+3i rmfs 3*pi/4 5*pi/4]);
s(88)=segment(3*M,[3+3i rmfs 5*pi/4 7*pi/4]);
s(89)=segment(3*M,[3+3i rmfs 7*pi/4 9*pi/4]);




% Now define all domains


d=domain.empty(29,0);

d(1)=domain(s([1 5 13 4]),[1 1 -1 -1]);
d(2)=domain(s([2 6 10 12 5]),[1 1 -1 -1 -1]);
d(3)=domain(s([3 7 9 6]),[1 1 -1 -1]);
d(4)=domain(s([13 12 11 16 14]),[1 1 1 -1 -1]);
d(5)=domain(s([10 9 8 18 11]),[1 1 1 -1 -1]);
d(6)=domain(s([16 17 23 15]),[1 1 -1 -1]);
d(7)=domain(s([18 19 21 17]),[1 1 -1 -1]);
d(8)=domain(s([23 22 26 24]),[1 1 -1 -1]);
d(9)=domain(s([21 20 27 22]),[1 1 -1 -1]);
d(10)=domain(s([26 30 35 41 29 25]),[1 1 -1 -1 1 -1]);
d(11)=domain(s([27 28 31 32 36 30]),[1 1 1 1 -1 -1]);
d(12)=domain(s([39 41 34 33]),[1 1 -1 -1]);
d(13)=domain(s([34 43 47 42]),[1 1 -1 -1]);
d(14)=domain(s([35 44 48 43]),[1 1 -1 -1]);
d(15)=domain(s([36 45 49 44]),[1 1 -1 -1]);
d(16)=domain(s([37 46 50 45]),[1 1 -1 -1]);
d(17)=domain(s([40 38 37 32]),[1 1 -1 -1]);
d(18)=domain(s([53 52 51 33 42 55 54]),[1 1 1 1 1 -1 -1]);
d(19)=domain(s([62 61 53 60]),[1 1 -1 -1]);
d(20)=domain(s([64 65 62 63]),[1 1 -1 -1]);
d(21)=domain(s([68 80 69 67 64 66]),[1 1 -1 1 -1 -1]);
d(22)=domain(s([70 71 75  73 72 80]),[1 1 -1 -1 -1 -1]);
d(23)=domain(s([75 76 78 74]),[1 1 -1 -1]);
d(24)=domain(s([78 79 57 77]),[1 1 -1 -1]);
d(25)=domain(s([57 58 59 46 38 81 56]),[1 1 -1 -1 -1 1 -1]);
d(26)=domain(s([55 47 48 49 50 59 82 86 83]),[1 1 1 1 1 1 1 1 -1]);
d(27)=domain(s([84 66 63 60 54 83 87]),[-1 1 1 1 1 1 1]);
d(28)=domain(s([88 85 70 68 84]),[1 -1 -1 -1 1]);
d(29)=domain(s([89 82 58 79 76 71 85]),[1 -1 -1 -1 -1 -1 1]);



ext=domain([],[],s([86 89 88 87]),-1);



% Make a list of all artificial segments and all boundary segments

bndindices=[1 2 3 4 14 15 24 25 7 8 19 20 28 ...
            29 39 51 52 61 65 67 69 72 73 74 77 81 40 31 56];
decompindices=setdiff(1:89,bndindices);
sbnd=s(bndindices);
sdecomp=s(decompindices);

% Now add basis functions to the domains

% First the singular corners

sopts.type='s'; sopts.cornermultipliers=[1 0 0 0]; sopts.rescale_rad=0;
d(1).addcornerbases(N,sopts);
sopts.cornermultipliers=[0 1 0 0 ];
d(3).addcornerbases(N,sopts);
sopts.cornermultipliers=[0 0 0 0 0 1];
d(10).addcornerbases(N,sopts);
sopts.cornermultipliers=[0 0 1 0 0 0];
d(11).addcornerbases(N,sopts);
sopts.cornermultipliers=[0 0 1 0 0 0 0];
d(18).addcornerbases(N,sopts);
sopts.cornermultipliers=[0 0 0 1 0 0];
d(21).addcornerbases(N,sopts);
sopts.cornermultipliers=[0 0 0 0 1 0];
d(22).addcornerbases(N,sopts);
sopts.cornermultipliers=[0 0 0 0 0 0 1];
d(25).addcornerbases(N,sopts);

% Now add regular Bessel functions to the rest

dsing=[1 3 10 11 18 21 22 25];
dnonsing=setdiff(1:29,dsing);

opts=struct('rescale_rad',0);
for j=1:length(dnonsing),
    dt=d(dnonsing(j));
    dt.addregfbbasis(dt.cloc(1),N,opts);
end

Z=@(t) 3+3i+3.2*sqrt(2)*exp(2i*pi*t);
Zp=@(t) 3.2*sqrt(2)*2i*pi*exp(2i*pi*t);
mfsopts=struct('fast',1,'nmultiplier',2,'eta',k);
ext.addmfsbasis({Z,Zp},N,mfsopts);

% Add the boundary conditions.

sdecomp.setmatch([k -k],[1 -1]);

% Need two cases depending on directions of segment

mpseglist=[65 61 39  2  8 19 20 40 81];
pseglist=[24 15 14 56 77 74];

s(mpseglist).setbc(-1,'D');
s(pseglist).setbc(1,'D');

% Define the problem

pr=scattering(ext,d);
pr.setoverallwavenumber(k);
pr.setincidentwave(-pi/2);

size(pr.A)
tic; pr.solvecoeffs; fprintf('\tcoeffs done in %.2g sec\n', toc)
fprintf('\tL2 bdry error norm = %g, coeff norm = %g\n', ...
    pr.bcresidualnorm, norm(pr.co))
resvec=[resvec,pr.bcresidualnorm];
    
    
o.bb=[-5 11 -5 11];
o.dx=0.05;

[ui gx gy] = pr.gridincidentwave(o);
u = pr.gridsolution(o);

figure;
imagesc(gx, gy, real(ui+u)); title('Full Field (Real Part)');
c = caxis; caxis([-1 1]*max(c));
axis equal tight;
colorbar;
set(gca,'ydir','normal'); 









