function avg = volumeAvg(field, V, mask)
%VOLUMEAVG  Volume-weighted average of FIELD over cells where MASK is
%true (cell volumes -- really per-unit-length areas -- vary with radius
%in this cylindrical grid, so a plain mean would be wrong).

avg = sum(field(mask).*V(mask)) / sum(V(mask));

end
