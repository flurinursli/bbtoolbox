
% get files folders
FolderInfo = dir('estimated_rupture*.txt')

nplanes = size(FolderInfo,1)

figure, hold on;

for ip = 1:nplanes
   
   FolderInfo(ip)
   %cell = importfile(FolderInfo(ip).name,1);
   cell = load(FolderInfo(ip).name);

   cell = cell';

   X = cell(2:3:end-1,:);
   Y = cell(1:3:end-1,:);
   Z = cell(3:3:end-1,:);
   C = cell(end,:);
  
   patch(X,Y,Z,C)

end

set(gca,'ZDir','Reverse')
grid on;
axis equal
   

FolderInfo = dir('connection*.txt')

nplanes = size(FolderInfo,1)

for ip = 1:nplanes

   FolderInfo(ip);
   %cell = importfile(FolderInfo(ip).name,2);
   cell = load(FolderInfo(ip).name);

   cell = cell';

   n = size(cell);

   for ic = 1:n(2)/2
     X = cell(2:3:12,ic);
     Y = cell(1:3:12,ic);
     Z = cell(3:3:12,ic);

     patch(X,Y,Z,'red')

     X = cell(2:3:12,n(2)/2 + ic);
     Y = cell(1:3:12,n(2)/2 + ic);
     Z = cell(3:3:12,n(2)/2 + ic);

     patch(X,Y,Z,'green')
   end

end

set(gca,'ZDir','Reverse')
grid on;
axis equal
