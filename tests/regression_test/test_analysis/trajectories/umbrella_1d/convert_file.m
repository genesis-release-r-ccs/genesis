center = 0:3:180;

for i = 1:numel(center)
  filename1 = sprintf('old/run_%d.dat', center(i));
  filename2 = sprintf('%d.dat', i);
  command = sprintf('cp %s %s', filename1, filename2)
  system(command);
    
  filename1 = sprintf('old/run_%d.out', center(i));
  filename2 = sprintf('%d.out', i);
  command = sprintf('cp %s %s', filename1, filename2)
  system(command);

  filename1 = sprintf('old/run_%d.nc', center(i));
  filename2 = sprintf('%d.dcd', i);
  [trj, box] = readnetcdf(filename1);
  writedcd(filename2, trj, box);
end

