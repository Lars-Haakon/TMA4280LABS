fileID = fopen('plotdata.txt','r');
m=fscanf(fileID,'%d', 1);
X=fscanf(fileID, '%f', m)';
Z=fscanf(fileID, '%f', [m m])';

figure;
mesh(X, X, Z);