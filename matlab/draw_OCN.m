function draw_OCN(OCN,X,varargin)

    % Create an inputParser object
    p = inputParser;

    % Define named properties and their default values
    addRequired(p, 'OCN');
    addRequired(p, 'X');
    addParameter(p, 'Draw_Borders', true);
    addParameter(p, 'Draw_River', true);
    addParameter(p, 'Draw_Frame', true);
    addParameter(p, 'Borders_Color', 'white');
    addParameter(p, 'River_Color', '#2776ea');
    addParameter(p, 'cmap', [0.03137254901960784 0.2549019607843137 0.3607843137254902; ...
            0.8 0.1607843137254902 0.21176470588235294]);
    addParameter(p, 'binary', false);

    % Parse the input arguments
    p.KeepUnmatched = true;
    parse(p,OCN,X,varargin{:})
    
    % Check if input is in grid format. If not, convert it.
    if all(size(X)==[OCN.nNodes,1])
        X = assignSC(OCN.MSC,X);
    end

    if all(size(X)==size(OCN.FD.A))
        Xtemp = X;
        X = zeros(size(OCN.MSC));
        for i = 1:length(Xtemp)
            X(OCN.geometry.XI(i),OCN.geometry.YI(i)) = Xtemp(i); 
        end
        clear Xtemp;
    end

    if p.Results.binary == false
        if ~isnan(X)
            h = imagesc(X','CDataMapping','direct','AlphaData',0.4);
            set(h, 'AlphaData', ~isnan(X'))
        end
    else
        imagesc(X','CDataMapping','direct','AlphaData',0.4);
        if all(X==0,'all')
            colormap([0.03137254901960784 0.2549019607843137 0.3607843137254902]);
        else
            colormap(p.Results.cmap)
        end
    end
    set(gca,'YDir','normal')
    hold on
    
    if p.Results.Draw_Borders == true
        drawborders(OCN.MSC,p.Results.Borders_Color)
    end

    if p.Results.Draw_River == true
        for i = 1:length(OCN.FD.A)
            if i ~= OCN.outlet
                if (OCN.FD.A(i)>=OCN.thrA)
                    line([OCN.FD.X(i)/OCN.cellsize OCN.FD.X(OCN.FD.downNode(i))/OCN.cellsize],...
                        [OCN.FD.Y(i)/OCN.cellsize OCN.FD.Y(OCN.FD.downNode(i))/OCN.cellsize],...
                            'color',p.Results.River_Color,...
                            'linewidth',0.9*(OCN.FD.A(i)*1e-1*1e-8)^0.3,...
                            'AlignVertexCenters','on')
                end
            end
        end
    end

    if p.Results.Draw_Frame == true
        line([0.5 0.5],[0.5 size(OCN.MSC,2)+.5],'color','k')
        line([size(OCN.MSC,1)+.5 size(OCN.MSC,1)+.5],[0.5 size(OCN.MSC,2)+.5],'color','k')
        line([0.5 size(OCN.MSC,1)+.5],[0.5 0.5],'color','k')
        line([0.5 size(OCN.MSC,1)+.5],[size(OCN.MSC,2)+.5 size(OCN.MSC,2)+.5],'color','k')

    end
    axis equal
    axis off
end
