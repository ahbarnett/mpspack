function addqprayleighbasis(d, seg, pm, varargin) %.......................
      % ADDQPRAYLEIGHBASIS - add Rayleigh basis to a domain, eg half-strip
      %
      % ADDREGQPRAYLEIGHBASIS(d, seg, pm, N) adds Rayleigh basis with order N,
      %  to domain d, radiative into the sense of the segment seg given by
      %  pm (=+-1). 
      %  For a qpstrip, if seg=pm=[], this is automatically into the
      %  half-strip's infinite direction.
      if isempty(seg), seg = d.seg; pm = d.pm; end % get T,B from half-strip 
      d.bas = {d.bas{:}, qprayleighbasis(seg, pm, varargin{:})};
      d.bas{end}.doms = d;    % tell basis it affects this domain
      %if isprop(d,'a')  % **** fails even when a is a property of domain!
        d.bas{end}.a = d.a;   % inherit Bloch phase (domain should have one)
      %else, warning('weird that domain does not have a Bloch phase'); end
  d.setupbasisdofs;  % ok since d will be some kind of unitcell
end
    
