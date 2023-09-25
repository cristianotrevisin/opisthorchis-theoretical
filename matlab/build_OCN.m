function OCN = build_OCN(name)

    load(fullfile('..','dataOCN', name));
    
    % Properties for drawing
    OCN.cellsize = 10000;%cellsize;
    cellsize = 10000;
    OCN.geometry.NX = ceil(max(X/cellsize)); 
    OCN.geometry.NY = ceil(max(Y/cellsize));

    
    
    % CCR = randi(25,[1 max(CTC)]);
    % 
    % CTD = zeros(1,length(CTC));
    % for i = 1:length(CTD)
    %     CTD(i) = CCR(CTC(i));
    % end
    
    %CTM = zeros(length(X),length(Y)); 
    
    OCN.MSC = zeros(OCN.geometry.NX,OCN.geometry.NY); % Map SubCatchment
    OCN.geometry.XI = ceil(X/cellsize);
    OCN.geometry.YI = ceil(Y/cellsize);

    OCN.geometry.SCX = SCX;OCN.geometry.SCY = SCY;
    for i = 1:length(X)
        % CTM(OCN.geometry.XI(i),OCN.geometry.YI(i)) = CTD(i);
        OCN.MSC(OCN.geometry.XI(i),OCN.geometry.YI(i)) = CTC(i); % CTC is pixel-to-catchment
    end

    OCN.nNodes = max(CTC,[],'all'); %numbers of nodes
    OCN.outlet = outlet;

    % To draw river networks
    OCN.FD.A = FD_A; OCN.FD.X = FD_X; OCN.FD.Y = FD_Y; OCN.FD.downNode = FD_downNode;

    % Attribute accumulating area for catchment -> assign to pixel
    OCN.SC_AccArea = A;
    % Area of subcatchments

    OCN.SC_Area = zeros(OCN.nNodes,1);
    for nn = 1:OCN.nNodes
        OCN.SC_Area(nn) = sum(CTC == nn)*cellsize^2;
    end


    % Build hydrological connectivity matrix

    OCN.W = zeros(SC.nNodes,SC.nNodes);

    for nn = 1:OCN.nNodes
        temp = find(downNode==nn);
        OCN.W(nn,temp) = 1;
    end


    % Get distances between nodes
    OCN.Dist = zeros(OCN.nNodes,OCN.nNodes);

    for sc1 = 1:OCN.nNodes
        for sc2 = 1:OCN.nNodes
            OCN.Dist(sc1,sc2) = sqrt((SCX(sc1)-SCX(sc2))^2+(SCY(sc1)-SCY(sc2))^2);
        end
    end

    % Rice fields suitability
    OCN.FD.Z = FD_Z;
    OCN.RicePaddy = rand(size(FD_Z)) <=  0.8*exp(-FD_Z*0.0008);
    OCN.SC_RicePaddy_Area = zeros(OCN.nNodes,1);
    for nn = 1:OCN.nNodes
        OCN.SC_RicePaddy_Area(nn) = sum(OCN.RicePaddy(CTC==nn))*cellsize^2;
    end

    % River specs
    OCN.SC_RiverLength = R_length;
    OCN.SC_RiverWidth = R_width;
    OCN.SC_RiverDepth = R_depth;

    OCN.SC_Volume = R_length.*R_width.*R_depth;

    OCN.FD.width = FD_width;
    OCN.FD.depth = FD_depth;

    % Connectivity
    OCN.distW = DSC+DSC'+ DSU+DSU';


end