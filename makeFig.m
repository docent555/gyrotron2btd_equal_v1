function [lHandleFmax, lHandleFabs, lHandleJmax, lHandleJabs, hF] = makeFig(zf, zp, t)

hF = figure;

axFmax = subplot(2,2,1);
axFmax.XLim = [0 t(end)];
lHandleFmax = line(axFmax, nan, nan); % Generate a blank line and return the line handle
lHandleFmax.Color = 'black';
lHandleFmax.LineWidth = 1;
axFmax.FontSize = 12;
axFmax.XLabel.String = 't';
axFmax.YLabel.String = '\Omega';
grid(axFmax,'on')

axFabs = subplot(2,2,3);
axFabs.XLim = [zf(1) zf(end)];
lHandleFabs = line(axFabs, nan, nan); % Generate a blank line and return the line handle
lHandleFabs.Color = 'black';
lHandleFabs.LineWidth = 1;
axFabs.FontSize = 12;
axFabs.XLabel.String = 'z';
axFabs.YLabel.String = '|F|';
lHandleFabs.XData = zf;
grid(axFabs,'on')

axJmax = subplot(2,2,2);
axJmax.XLim = [0 t(end)];
lHandleJmax = line(axJmax, nan, nan); % Generate a blank line and return the line handle
lHandleJmax.Color = 'black';
lHandleJmax.LineWidth = 1;
axJmax.FontSize = 12;
axJmax.XLabel.String = 't';
axJmax.YLabel.String = 'Efficiency';
grid(axJmax,'on')

axJabs = subplot(2,2,4);
axJabs.XLim = [zp(1) zp(end)];
lHandleJabs = line(axJabs, nan, nan); % Generate a blank line and return the line handle
lHandleJabs.Color = 'black';
lHandleJabs.LineWidth = 1;
axJabs.FontSize = 12;
axJabs.XLabel.String = 'z';
axJabs.YLabel.String = '|J|';
lHandleJabs.XData = zp;
grid(axJabs,'on')
end