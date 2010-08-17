% WEYLCOUNTCHECK - basic Weyl law check for missing Dirichlet eigenwavenumbers
%
% function [k_weyl] = weylchoutcheck(k_lo, ks, perim, area {, dkfig {, j_lo}})
%
% test (new figure) array of wavenumbers against the Weyl law (first 2 terms).
% k_lo is where the sequence started. perim term for Dirichlet BCs.
% k_weyl is sequence of Weyl-computed k values, assuming lowest is j_lo
%
% Copied from old Weyl code from Courant days.
% barnett 11/24/03
% 2/1/04 included figure for talk: dkfig causes to show N(k), N_weyl(k).
% 4/12/06 added k_weyl

function [k_weyl] = weylcountcheck(k_lo, k, perim, area, dkfig)

i = reshape(1:length(k), size(k));
a = area/(4*pi);
b = -perim/(4*pi);

%figure;
plot(k, a*k.*k + b*k - i+0.5);
avg = mean(a*k.*k + b*k - i+0.5);
%hold on; hline(avg-2); hline(avg-1); hline(avg); hline(avg+1); hline(avg+2);
hold off;
%title('Weyl test of k sequence (jump up = missed state)');
set(gca, 'fontsize', 14);
xlabel('k'); ylabel('N_{Weyl}(k) - N(k)');

if nargin>4 && dkfig>0
  figure;
  ks = k_lo:dkfig/200:k_lo+dkfig;
  for i=1:length(ks)
    kf = ks(i);
    Nk(i) = length(find(k<kf));
    Nwk(i) = a*kf*kf + b*kf - (a*k_lo*k_lo + b*k_lo);
  end
  plot(ks, Nk, '-',  ks, Nwk, '--');
  xlabel('k'); ylabel('N');
  axis([k_lo k_lo+dkfig 0 Nk(length(ks))]);
  %legend('N(k)', 'N_{Weyl}(k)');

end

if nargin<5
  j_lo = 1;
end

% solve quadratic for k(N):
k_weyl = (-b +sqrt(b^2 + 4*a*((1:numel(k))+j_lo-1)))/(2*a);
