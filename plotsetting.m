function plotsetting(hfig)
%% In most situations, the lines below are default. However, you can play
%% around with some of the parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 1; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),...
    'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',...
    [3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','manual','PaperUnits','centimeters',...
    'PaperSize',[pos(3), pos(4)])
set(hfig, 'PaperPosition',[0 0 pos(3) pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
