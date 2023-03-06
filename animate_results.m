function animate_results(Y,MSC,FD,thrA,cellsize)


% for t = 1:total_time
%     t
%     figure(4)
%     imshow(assignSC(MSC,WH(t,:))')
%     set(gca,'YDir','normal')
%     colorbar
% 
% end

% h = figure;
% 
% axis tight manual
% ax = gca;
% ax.NextPlot = 'replaceChildren';
% 
% 
% M(total_time) = struct('cdata',[],'colormap',[]);
% 
% h.Visible = 'off';
% 
% for t = 1:total_time
% 
%     subplot(2,2,1)
% 
%     image(assignSC(MSC,WH(t,:))')
%     set(gca,'YDir','normal')
%     title('Worm burden in humans')
% 
%     subplot(2,2,2)
% 
%     image(assignSC(MSC,SI(t,:))')
%     set(gca,'YDir','normal')
%     title('Egg burden in snails')
% 
%     subplot(2,2,3)
% 
%     image(assignSC(MSC,FI(t,:))')
%     set(gca,'YDir','normal')
%     title('Cyst burden in fish')
% 
%     subplot(2,2,4)
%     title(['t = ',num2str(t)])
% 
%     drawnow
% 
%     M(t) = getframe;
% end
% 
% h.Visible = 'on';
% 
% movie(M);

%  h = figure;
% 
% axis tight manual
% ax = gca;
% ax.NextPlot = 'replaceChildren';
% 
% M(total_time) = struct('cdata',[],'colormap',[]);
% 
% h.Visible = 'off';
% 
% for t = 1:total_time
% 
%     image(assignSC(MSC,WH(t,:))')
%     set(gca,'YDir','normal')
%     title(t)
% 
%     drawnow
% 
%     M(t) = getframe;
% end
% 
% h.Visible = 'on';
% 
% movie(M);
h = figure;

for t = 1:50:size(Y,1)
t
    imagesc(assignSC(MSC,Y(t,:)))
    set(gca,'YDir','normal')
    %drawborders(MSC)
    colorbar
    for i = 1:FD.nNodes
        if i ~= FD.outlet
            if (FD.A(i)>=thrA && ...
                    abs(FD.X(i)-FD.X(FD.downNode(i)))<=cellsize && ...
                    abs(FD.Y(i)-FD.Y(FD.downNode(i)))<=cellsize)
                line([FD.X(i) FD.X(FD.downNode(i))],...
                    [FD.Y(i) FD.Y(FD.downNode(i))],...
                    'linewidth', 0.5+4.5*sqrt(FD.A(i)/(FD.nNodes*cellsize^2)),...
                    'color','white')
            end
        end
    end
    title(t)

    exportgraphics(gcf,'test.gif','Append',true);
end


function M = assignSC(MSC,X)
    M = zeros(size(MSC));
    for iz = 1:length(X)
        M(MSC==iz) = X(iz);
    end
end


end