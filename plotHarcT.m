function plotHarcT(t, j, x, color)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlab M-file                Author: Yuchun Li
%
% Description: yet another function for ploting hybrid arcs projected in the t-axis
% parameter dSample = 200 specifies the rate of down sampling, if there are more than
% dSample points for a line segment, it will resample the points according
% to the dSample rate, in turn, it generates a plot with much smaller file
% size for many cases. 
% Note: it will not work if the reason for too many data points is because
% of high frequency jumps. 
%
% sample usage: plotHarcT(t, j, x, 'r');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% max (t,j) in simulation
T = t(end);
J  = j(end);
% dSample = 200;
dSample = 20;

for i = 0:J
    hold on
    indexi = find(j == i);
    if (length(indexi) > dSample)
        step = round(length(indexi) / dSample);
        indexiNew = indexi(1:step:end);
        if indexiNew(end) ~= indexi(end)
            indexiNew(end+1) = indexi(end);
        end
    else
        indexiNew = indexi;
    end
    tPost = t(indexiNew(1));
    xPost = x(indexiNew(1));
    
    if (i > 0) 
        plot([tPre tPost], [xPre xPost], 'r--');
    end
    
    plot(t(indexiNew),x(indexiNew),color, 'linewidth', 2);
    tPre = t(indexiNew(end));
    xPre = x(indexiNew(end));
end

xlabel('t')
grid on
set(gca,'Box','on')
set(gca,'FontSize',10)