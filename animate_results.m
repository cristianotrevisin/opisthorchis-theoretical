function animate_results(y,MSC)


total_time = size(y,1);

WH = y(:,1:3:end);

SI = y(:,2:3:end);

FI = y(:,3:3:end);

for t = 1:total_time
    t
    imshow(assignSC(MSC,WH(t,:))')

end

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


function M = assignSC(MSC,X)
    M = zeros(size(MSC));
    for iz = 1:length(X)
        M(MSC==iz) = X(iz);
    end
end


end