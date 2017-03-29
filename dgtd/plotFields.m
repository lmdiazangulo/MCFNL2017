function plotFields(e, time)
% function plotFields(e)

% --------- Preliminar computations ---------------------------------------
constants;
K = numel(e);
Np = size(e(1).x,1);
N = Np - 1;
% --------- Plots Electric field ------------------------------------------
subplot(2,1,1);
for k = 1:K
    plot( e(k).x , e(k).E ); hold on;
    plot([ e(k).x(1) , e(k).x(1) ],[-1.1, 1.1],'k');
%     set(gcf, 'Color', [1 0 0]);
end
grid on;
axis([e(1).x(1) e(K).x(Np) -1.1 +1.1]);
ylabel('E field [V/m]');
xlabel('Distance in meters');
title(sprintf('DGTD-P^%d, time = %.2f ns',N, time*1e9), ...
              'Color',[0 0 1],'FontSize', 10)
hold off;
% --------- Plots Magnetic field ------------------------------------------
subplot(2,1,2);
for k = 1:K
    plot( e(k).x, e(k).H ); hold on;
    plot([ e(k).x(1), e(k).x(1) ],[-1.1, 1.1], 'k');
end
grid on;
axis([e(1).x(1) e(K).x(Np) -1.1/e(1).Z +1.1/e(1).Z]);
ylabel('H field [A/m]');
xlabel('Distance in meters');
hold off;
% -------- Shows as a motion picture fashion ------------------------------
drawnow;
