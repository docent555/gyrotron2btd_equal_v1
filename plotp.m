function plotp(x,p)
figure;
% plot(x,abs(p(:,end)))
for i = 1:size(p,2)
    plot(x,abs(p(:,i))); hold on
end
hold off
end

