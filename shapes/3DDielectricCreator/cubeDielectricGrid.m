
function [eps_final, interiorCoords, borderCoords] = cubeDielectricGrid(cellx, celly,...
    cellz, SingleCellDims,...
    inclusionSize, diel)
    %SingleCellDims has [dimx, dimy, dimz]
    
    Lx = SingleCellDims(1);
    Ly = SingleCellDims(2);
    Lz = SingleCellDims(3);
    
    numCells = cellx*celly*cellz;
    %% extract cubeface coordinates
    cubeface = cell(3,1);
    xy = combvec(2:Lx, 2:Ly).'; %% all possible permutations
    xz = combvec(1:Lx, 1:Lz).';
    yz = combvec(2:Ly, 1:Lz).';

    cubeface{1} = [xy, ones(length(xy),1)];
    cubeface{2} = [ones(length(yz),1), yz];
    cubeface{3} = [xz(:,1), ones(length(xz),1), xz(:,2)];

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

          end
       end
    end

    %% now we need to wrap the edges
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

    

end