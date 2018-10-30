Lx = 10; Ly = Lx; Lz = Lx;
SingleCellDims= [Lx, Ly, Lz];
diel = 12;
inclusionSize = 7;
cellx = 2; celly = 2; cellz = 2;
numCells = cellx*celly*cellz;
%% extract cubeface coordinates
cubeface = cell(3,1);
xy = combvec(2:Lx, 2:Ly).';
xz = combvec(1:Lx, 1:Lz).';
yz = combvec(2:Ly, 1:Lz).';

cubeface{1} = [xy, ones(length(xy),1)];
cubeface{2} = [ones(length(yz),1), yz];
cubeface{3} = [xz(:,1), ones(length(xz),1), xz(:,2)];

s = 60;
for i = 1:3
   a = cubeface{i};
   scatter3(a(:,1), a(:,2), a(:,3), 'filled')
   hold on;
   s = s-20;
end


numFC = length(xz)+length(yz)+length(xy);


%% interior coords
xyz = combvec(2:Lx, 2:Ly, 2:Lz).';

eps_cells = cell(cellx, celly, cellz);
borderCoords = cell(cellx, celly, cellz);
interiorCoords = cell(cellx, celly, cellz);
prevcoord = [1,1,1];
for i = 1:cellx
   for j = 1:celly
      for k = 1:cellz
         %% face coordinates or border
         facecoords = cell2mat(cubeface);
         interiorcoords = xyz;
         xdisp = (i-1)*Lx; ydisp = (j-1)*Ly; zdisp = (k-1)*Lz;
         disp = [xdisp*ones(numFC,1), ydisp*ones(numFC,1), zdisp*ones(numFC,1)];
         facecoords = facecoords+disp;
         
         %% create hte interoir
         eps_proto = ones(SingleCellDims);
         eps_proto = createCube(eps_proto, inclusionSize, diel);
         eps_cells{i,j,k} = eps_proto;
         borderCoords{i,j,k} = facecoords;
         %% interior Coords
         numIC = length(xyz);
         interiordisp = [xdisp*ones(numIC,1), ydisp*ones(numIC,1), zdisp*ones(numIC,1)];
         interiorcoords = xyz + interiordisp;
         interiorCoords{i,j,k} = interiorcoords
         
         %%diagnose
%          scatter3(interiorcoords(:,1), interiorcoords(:,2), interiorcoords(:,3), 'filled');
%          hold on;
%          scatter3(facecoords(:,1), facecoords(:,2), facecoords(:,3), 'filled');
%          hold on;
      end
   end
end

%% now we need to wrap the edges
s = 60 ;
facewrap = cell(3,1);
Ltx = cellx*Lx+1; Lty = celly*Ly+1; Ltz = cellz*Lz+1;
xy = combvec(1:Ltx-1, 1:Lty-1).';
xz = combvec(1:Ltx, 1:Ltz).';
yz = combvec(1:Lty-1, 1:Ltz).';

facewrap{1} = [xy, Ltz*ones(length(xy),1)];
facewrap{2} = [Ltx*ones(length(yz),1), yz];
facewrap{3} = [xz(:,1), Lty*ones(length(xz),1), xz(:,2)];
borderCoords = reshape(borderCoords, [numCells, 1])

for i = 1:3
   borderCoords{numCells+i} = facewrap{i}; 
end
eps_final = cell2mat(eps_cells);
eps_final(Ltx,1,1) = 1;
eps_final(1,Lty,1) = 1;
eps_final(1,1,Ltz) = 1;
eps_final(eps_final == 0) = 1;

for i = 1:3
   a = facewrap{i};
   scatter3(a(:,1), a(:,2), a(:,3), 'filled')
   hold on;
   s = s-20;
end

[eps_final, interiorCoords, borderCoords] = cubeDielectricGrid(cellx, celly,...
    cellz, SingleCellDims,...
    inclusionSize, diel)

borders = vertcat(borderCoords{:});
interiors = vertcat(interiorCoords{:});

figure;
scatter3(borders(:,1), borders(:,2), borders(:,3), 'filled')
figure; 
scatter3(interiors(:,1), interiors(:,2), interiors(:,3), 'filled')

%% DUAL SEPARATOR
% Lx = 10; Ly = Lx; Lz = Lx;
% SingleCellDims= [Lx, Ly, Lz];
% diel = 12;
% inclusionSize = 2;
% cellx = 2; celly = 2; cellz = 2;
% 
% %% extract cubeface coordinates
% cubeface = cell(6,1);
% xy = combvec(1:Lx, 1:Ly).'
% numFC = 6*length(xy);
% b1 = ones(numFC/6, 1);
% cubeface{3} = [xy, b1]; cubeface{4} = [xy, Lx*b1];
% cubeface{1} = [b1, xy]; cubeface{5} = [Lx*b1, xy];
% cubeface{2} = [xy(:,1), b1, xy(:,2)];
% cubeface{6} = [xy(:,1), Lx*b1, xy(:,2)];

%%% interior coords
% xyz = combvec(2:Lx-1, 2:Ly-1, 2:Lz-1).';
% 
% eps_cells = cell(cellx, celly, cellz);
% borderCoords = cell(cellx, celly, cellz);
% interiorCoords = cell(cellx, celly, cellz);
% prevcoord = [1,1,1];
% for i = 1:cellx
%    for j = 1:celly
%       for k = 1:cellz
%          %% face coordinates or border
%          facecoords = cell2mat(cubeface);
%          interiorcoords = xyz;
%          xdisp = (i-1)*Lx; ydisp = (j-1)*Ly; zdisp = (k-1)*Lz;
%          disp = [xdisp*ones(numFC,1), ydisp*ones(numFC,1), zdisp*ones(numFC,1)];
%          facecoords = facecoords+disp;
%          
%          %% create hte interoir
%          eps_proto = ones(SingleCellDims);
%          eps_proto = createCube(eps_proto, inclusionSize, diel);
%          eps_cells{i,j,k} = eps_proto;
%          borderCoords{i,j,k} = facecoords;
%          %% interior Coords
%          numIC = length(xyz);
%          interiordisp = [xdisp*ones(numIC,1), ydisp*ones(numIC,1), zdisp*ones(numIC,1)];
%          interiorcoords = xyz + interiordisp;
%          interiorCoords{i,j,k} = interiorcoords
%          
%          %%diagnose
%          scatter3(interiorcoords(:,1), interiorcoords(:,2), interiorcoords(:,3), 'filled');
%          hold on;
%          scatter3(facecoords(:,1), facecoords(:,2), facecoords(:,3), 'filled');
%          hold on;
%       end
%    end
% end
% eps_final = cell2mat(eps_cells);
% %% generate the cube face coordinates


