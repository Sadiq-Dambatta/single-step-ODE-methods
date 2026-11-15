%% Haruna Sadiq Sadiq - Comparative Error and Convergence Analysis
% FULL file: Tables 1,2,3,4 + Sec 6.3 + YOUR Figures 1-4 code + YOUR helpers
% Run: main_complete.m

clear; clc; close all;
P = get_problems();
h_list = [0.2, 0.1, 0.05, 0.025, 0.0125];
allMaxErr = zeros(numel(P), 4, length(h_list));

% ============ PART A: TABLES 3 + TABLE 4 (Fitted orders) ============
for k = 1:numel(P)
    for j = 1:length(h_list)
        h = h_list(j);
        [xe, ye] = euler(P(k).f, P(k).x0, P(k).y0, h, P(k).xend);
        [~, ym] = modified_euler(P(k).f, P(k).x0, P(k).y0, h, P(k).xend);
        [~, yt] = taylor2(P(k).f, P(k).df, P(k).x0, P(k).y0, h, P(k).xend);
        [~, yr] = rk4(P(k).f, P(k).x0, P(k).y0, h, P(k).xend);
        yex = P(k).exact(xe);
        allMaxErr(k,1,j)=max(abs(yex-ye)); allMaxErr(k,2,j)=max(abs(yex-ym));
        allMaxErr(k,3,j)=max(abs(yex-yt)); allMaxErr(k,4,j)=max(abs(yex-yr));
    end
    fprintf('P%d Fitted orders: ', k);
    for m = 1:4
        p_fit = polyfit(log(h_list), log(squeeze(allMaxErr(k,m,:))'), 1);
        fprintf('%.2f ', p_fit(1));
    end
    fprintf('\n');
end

% ============ PART B: TABLE 2 - VARIATIONAL EQUATION (Sec 5) ============
fprintf('\n========== TABLE 2 ==========\n');
fprintf('%8s %10s %10s %10s %10s\n','Problem','C_Heun','C_Taylor','PredRatio','ObsRatio');
D = get_derivs(); he = 1e-4; te = 0:he:1;
for k = 1:numel(P)
    yex = P(k).exact(te); fval = P(k).f(te, yex);
    ft = D(k).ft(te,yex); fy = D(k).fy(te,yex);
    ftt = D(k).ftt(te,yex); fty = D(k).fty(te,yex); fyy = D(k).fyy(te,yex);
    yppp = ftt + 2*fval.*fty + fval.^2.*fyy + ft.*fy + fval.*fy.^2;
    cT = yppp/6;
    cH = (-(ftt + 2*fval.*fty + fval.^2.*fyy) + 2*fy.*(ft + fval.*fy))/12;
    eH = solve_variational(te, fy, cH, he); eT = solve_variational(te, fy, cT, he);
    C_H = max(abs(eH)); C_T = max(abs(eT));
    fprintf('%8d %10.4f %10.4f %10.2f %10.2f\n', k, C_H, C_T, C_T/C_H, allMaxErr(k,3,5)/allMaxErr(k,2,5));
end

% ============ PART C: SECTION 6.3 - COST ============
fprintf('\n========== SECTION 6.3 ==========\n');
target = allMaxErr(1,4,1); fprintf('Target = %.3e\n', target);
methods = {'Euler','Heun','Taylor2','RK4'};
for m = 1:4
    N = find_min_N(m, P(1), target);
    if m==1, feval=N*1; elseif m==2, feval=N*2; elseif m==3, feval=N*1; else, feval=N*4; end
    fprintf('%s: N=%d, f-evals=%d\n', methods{m}, N, feval);
end

% ============ PART D: YOUR FIGURES 1-4 CODE (VERBATIM) ============
%% Figures 1-4
h_list = [0.2, 0.1, 0.05, 0.025, 0.0125];
P      = get_problems();
Emax   = zeros(numel(P), numel(h_list), 4);
for k = 1:numel(P)
    for j = 1:numel(h_list)
        h = h_list(j);
        x0 = P(k).x0;  y0 = P(k).y0;  xend = P(k).xend;
        [xe, ye] = euler(P(k).f, x0, y0, h, xend);
        [~,  ym] = modified_euler(P(k).f, x0, y0, h, xend);
        [~,  yt] = taylor2(P(k).f, P(k).df, x0, y0, h, xend);
        [~,  yr] = rk4(P(k).f, x0, y0, h, xend);
        yex = P(k).exact(xe);
        E = [abs(yex-ye); abs(yex-ym); abs(yex-yt); abs(yex-yr)];
        Emax(k,j,:) = max(E, [], 2);
    end
end

SAVE_FIGS = true;
all4 = 1:4;
refs = [1 1; 2 2; 4 4];   % reference slopes: [order, method index]
conv_figure(h_list, Emax, 1, all4,  refs,       'fig1', SAVE_FIGS);  % Fig. 1
conv_figure(h_list, Emax, 2, [2 3], zeros(0,2), 'fig2', SAVE_FIGS);  % Fig. 2
solution_figure(P(2), 0.2, 'fig3', SAVE_FIGS);                       % Fig. 3
conv_figure(h_list, Emax, 7, all4,  refs,       'fig4', SAVE_FIGS);  % Fig. 4

% ============ FUNCTIONS ============
function [x, y] = euler(f, x0, y0, h, xend)
    N = round((xend - x0)/h); x = x0 + (0:N)*h; y = zeros(1,N+1); y(1)=y0;
    for n=1:N, y(n+1)=y(n)+h*f(x(n),y(n)); end
end
function [x, y] = modified_euler(f, x0, y0, h, xend)
    N = round((xend - x0)/h); x = x0 + (0:N)*h; y = zeros(1,N+1); y(1)=y0;
    for n=1:N, k1=f(x(n),y(n)); yp=y(n)+h*k1; k2=f(x(n+1),yp); y(n+1)=y(n)+(h/2)*(k1+k2); end
end
function [x, y] = taylor2(f, df, x0, y0, h, xend)
    N = round((xend - x0)/h); x = x0 + (0:N)*h; y = zeros(1,N+1); y(1)=y0;
    for n=1:N, y(n+1)=y(n)+h*f(x(n),y(n))+(h^2/2)*df(x(n),y(n)); end
end
function [x, y] = rk4(f, x0, y0, h, xend)
    N = round((xend - x0)/h); x = x0 + (0:N)*h; y = zeros(1,N+1); y(1)=y0;
    for n=1:N
        k1=f(x(n),y(n)); k2=f(x(n)+h/2,y(n)+(h/2)*k1);
        k3=f(x(n)+h/2,y(n)+(h/2)*k2); k4=f(x(n)+h,y(n)+h*k3);
        y(n+1)=y(n)+(h/6)*(k1+2*k2+2*k3+k4);
    end
end
function e = solve_variational(t, fy, c, h)
    e = zeros(size(t));
    for n=1:length(t)-1, e(n+1)=e(n)+h*(fy(n)*e(n)+c(n)); end
end
function Nmin = find_min_N(method, prob, target)
    lo=1; hi=100000;
    while lo<hi
        mid=floor((lo+hi)/2); h=1/mid;
        if method==1, [~,y]=euler(prob.f,prob.x0,prob.y0,h,prob.xend);
        elseif method==2, [~,y]=modified_euler(prob.f,prob.x0,prob.y0,h,prob.xend);
        elseif method==3, [~,y]=taylor2(prob.f,prob.df,prob.x0,prob.y0,h,prob.xend);
        else, [~,y]=rk4(prob.f,prob.x0,prob.y0,h,prob.xend); end
        x = prob.x0 + (0:mid)*h; Emax=max(abs(prob.exact(x)-y));
        if Emax<=target, hi=mid; else, lo=mid+1; end
    end
    Nmin=lo;
end
function conv_figure(h, Emax, k, midx, refs, stem, doSave)
    labels = {'Euler','Modified Euler (Heun)','Taylor-2','RK4'};
    styles = {'-o','-s','--d','-^'};
    colors = [0.85 0.33 0.10; 0.00 0.45 0.74; 0.47 0.67 0.19; 0 0 0];
    fig = figure('Color','w','Units','centimeters','Position',[2 2 12 9]);
    ax  = axes(fig);  hold(ax,'on');
    for m = midx
        if m == 3, mfc = 'none'; else, mfc = colors(m,:); end
        plot(ax, h, squeeze(Emax(k,:,m)), styles{m}, 'Color', colors(m,:), ...
             'LineWidth', 1.3, 'MarkerSize', 5, 'MarkerFaceColor', mfc, ...
             'DisplayName', labels{m});
    end
    for r = 1:size(refs,1)
        p  = refs(r,1);  mi = refs(r,2);
        e0 = 0.5*Emax(k,end,mi);
        hh = [h(end) h(1)];
        plot(ax, hh, e0*(hh/h(end)).^p, 'k:', 'LineWidth', 0.8, 'HandleVisibility','off');
        text(ax, h(end)*1.05, e0*0.3, sprintf('O(h^{%d})', p), 'FontSize', 8);
    end
    set(ax, 'XScale','log', 'YScale','log');
    hs = sort(h);
    set(ax, 'XTick', hs, 'XTickLabel', arrayfun(@(v) sprintf('%g',v), hs, 'UniformOutput', false));
    xlim(ax, [min(h)*0.8, max(h)*1.25]);
    xlabel(ax, 'Step size h');  ylabel(ax, 'Maximum absolute error');
    grid(ax,'on');  box(ax,'on');
    legend(ax, 'Location', 'northwest');
    if doSave, save_figure(fig, stem); end
end
function solution_figure(prob, h, stem, doSave)
    [xe, ye] = euler(prob.f, prob.x0, prob.y0, h, prob.xend);
    [xr, yr] = rk4(prob.f, prob.x0, prob.y0, h, prob.xend);
    xx  = linspace(prob.x0, prob.xend, 200);
    fig = figure('Color','w','Units','centimeters','Position',[2 2 12 9]);
    ax  = axes(fig);  hold(ax,'on');
    plot(ax, xx, prob.exact(xx), 'k-', 'LineWidth', 1.3, 'DisplayName', 'Exact');
    plot(ax, xe, ye, 'ro', 'MarkerFaceColor','r', 'DisplayName', sprintf('Euler (h = %.1f)', h));
    plot(ax, xr, yr, 'bs', 'MarkerFaceColor','b', 'DisplayName', sprintf('RK4 (h = %.1f)', h));
    xlabel(ax, 't');  ylabel(ax, 'y(t)');
    grid(ax,'on');  box(ax,'on');
    legend(ax, 'Location', 'northwest');
    if doSave, save_figure(fig, stem); end
end
function save_figure(fig, stem)
    try
        exportgraphics(fig, [stem '.pdf'], 'ContentType', 'vector');
    catch
        pos = get(fig, 'Position');
        set(fig, 'PaperUnits','centimeters', 'PaperSize', pos(3:4), 'PaperPosition', [0 0 pos(3:4)]);
        print(fig, stem, '-dpdf', '-painters');
    end
end
function D = get_derivs()
    D(1).ft=@(t,y) 1+0*t; D(1).fy=@(t,y) 1+0*t; D(1).ftt=@(t,y) 0*t; D(1).fty=@(t,y) 0*t; D(1).fyy=@(t,y) 0*t;
    D(2).ft=@(t,y) y; D(2).fy=@(t,y) t; D(2).ftt=@(t,y) 0*t; D(2).fty=@(t,y) 1+0*t; D(2).fyy=@(t,y) 0*t;
    D(3).ft=@(t,y) -2*t; D(3).fy=@(t,y) 1+0*t; D(3).ftt=@(t,y) -2+0*t; D(3).fty=@(t,y) 0*t; D(3).fyy=@(t,y) 0*t;
    D(4).ft=@(t,y) -sin(t); D(4).fy=@(t,y) -1+0*t; D(4).ftt=@(t,y) -cos(t); D(4).fty=@(t,y) 0*t; D(4).fyy=@(t,y) 0*t;
    D(5).ft=@(t,y) 0*t; D(5).fy=@(t,y) -0.5+0*t; D(5).ftt=@(t,y) 0*t; D(5).fty=@(t,y) 0*t; D(5).fyy=@(t,y) 0*t;
    D(6).ft=@(t,y) 0*t; D(6).fy=@(t,y) 1-2*y; D(6).ftt=@(t,y) 0*t; D(6).fty=@(t,y) 0*t; D(6).fyy=@(t,y) -2+0*t;
    D(7).ft=@(t,y) -2*y.^2; D(7).fy=@(t,y) -4*t.*y; D(7).ftt=@(t,y) 0*t; D(7).fty=@(t,y) -4*y; D(7).fyy=@(t,y) -4*t;
end
function P = get_problems()
    P(1).name='P1: dy/dx = x + y,  y(0)=1'; P(1).f=@(x,y) x+y; P(1).df=@(x,y) 1+x+y; P(1).exact=@(x) 2*exp(x)-x-1; P(1).x0=0; P(1).y0=1; P(1).xend=1;
    P(2).name='P2: dy/dx = x*y,  y(0)=1'; P(2).f=@(x,y) x.*y; P(2).df=@(x,y) y.*(1+x.^2); P(2).exact=@(x) exp(x.^2/2); P(2).x0=0; P(2).y0=1; P(2).xend=1;
    P(3).name='P3: dy/dx = y - x^2 + 1,  y(0)=0.5'; P(3).f=@(x,y) y-x.^2+1; P(3).df=@(x,y) y-x.^2-2*x+1; P(3).exact=@(x) x.^2+2*x+1-0.5*exp(x); P(3).x0=0; P(3).y0=0.5; P(3).xend=1;
    P(4).name='P4: dy/dx = cos(x) - y,  y(0)=0'; P(4).f=@(x,y) cos(x)-y; P(4).df=@(x,y) -sin(x)-cos(x)+y; P(4).exact=@(x) 0.5*(sin(x)+cos(x)-exp(-x)); P(4).x0=0; P(4).y0=0; P(4).xend=1;
    P(5).name='P5: dT/dt = -0.5(T-25),  T(0)=100'; P(5).f=@(x,y) -0.5*(y-25); P(5).df=@(x,y) 0.25*(y-25); P(5).exact=@(x) 25+75*exp(-0.5*x); P(5).x0=0; P(5).y0=100; P(5).xend=1;
    P(6).name='P6: dy/dx = y(1-y),  y(0)=0.1'; P(6).f=@(x,y) y.*(1-y); P(6).df=@(x,y) (1-2*y).*y.*(1-y); P(6).exact=@(x) 1./(1+9*exp(-x)); P(6).x0=0; P(6).y0=0.1; P(6).xend=1;
    P(7).name='P7: dy/dx = -2xy^2,  y(0)=1'; P(7).f=@(x,y) -2*x.*y.^2; P(7).df=@(x,y) -2*y.^2+8*x.^2.*y.^3; P(7).exact=@(x) 1./(1+x.^2); P(7).x0=0; P(7).y0=1; P(7).xend=1;
end
