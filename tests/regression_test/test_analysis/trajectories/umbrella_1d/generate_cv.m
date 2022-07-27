index = [9 15 17 19];
center = 0:3:180;

for i = 1:numel(center)
  filename = sprintf('%d.dcd', i)
  trj = readdcd(filename);
  phi = calcdihedral(trj, index).*180./pi;

  filename = sprintf('%d.dat', i);
  fid = fopen(filename, 'w');
  for istep = 1:numel(phi)
    fprintf(fid, '%d %f\n', istep*5000, phi(istep));
  end
  fclose(fid);
end

