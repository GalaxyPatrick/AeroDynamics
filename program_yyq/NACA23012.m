%for plotting and discretizing NACA23012

%%%%%%%%%%%%%%%%%%%%%%%parameter setting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global m k1 t x gap_smallx gap_largex cpoint ...
       need_plot_continuous need_plot_discrete
m = 0.2025;
k1 = 15.957;
t = 12/100;

x = 0:1e-3:1; %total number of interval is 1000
cpoint = 0.04; %to seperate small x and large x
gap_smallx = 4;
gap_largex = 40;
need_plot_continuous = 0;
need_plot_discrete = 1;

%%%%%%%%%%%%%%%%%%%%%%%main function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

main_23012();

function main_23012()
    global m k1 t x gap_smallx gap_largex cpoint ...
           need_plot_continuous need_plot_discrete

    if need_plot_continuous
        x_u = arrayfun(@(x) get_x_u(x, m, k1, t), x);
        x_l = arrayfun(@(x) get_x_l(x, m, k1, t), x);
        y_u = arrayfun(@(x) get_y_u(x, m, k1, t), x);
        y_l = arrayfun(@(x) get_y_l(x, m, k1, t), x);
        plot_23012(x_u, x_l, y_u, y_l, 'c');
        hold on
    end
    if need_plot_discrete
        x_d = discretizer(x, cpoint, gap_smallx, gap_largex);
        x_u = arrayfun(@(x) get_x_u(x, m, k1, t), x_d);
        x_l = arrayfun(@(x) get_x_l(x, m, k1, t), x_d);
        y_u = arrayfun(@(x) get_y_u(x, m, k1, t), x_d);
        y_l = arrayfun(@(x) get_y_l(x, m, k1, t), x_d);
        plot_23012(x_u, x_l, y_u, y_l, 'd');
        hold on

        x_u_c = centerize(x_u);
        x_l_c = centerize(x_l);
        y_u_c = centerize(y_u);
        y_l_c = centerize(y_l);
        scatter(x_u_c, y_u_c, 18, [0.4660 0.6740 0.1880]);
        hold on
        scatter(x_l_c, y_l_c, 18, [0.4660 0.6740 0.1880]);

        writematrix(x_u, 'coord/x_u.csv');
        writematrix(x_u_c, 'coord/x_u_c.csv');
        writematrix(x_l, 'coord/x_l.csv');
        writematrix(x_l_c, 'coord/x_l_c.csv');
        writematrix(y_u, 'coord/y_u.csv');
        writematrix(y_u_c, 'coord/y_u_c.csv');
        writematrix(y_l, 'coord/y_l.csv');
        writematrix(y_l_c, 'coord/y_l_c.csv');
        writematrix(centerize(x_d), 'coord/x_d_c.csv');
    end
end

%%%%%%%%%%%%%%%%%plot and discretize function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_23012(x_u, x_l, y_u, y_l, c_or_d)
    if c_or_d == 'c'
        color = '#0072BD';
    else
        color = '#D95319';
    end
    title("NACA 23012");
    plot(x_u, y_u, 'Color', color, 'LineWidth', 1);
    hold on;
    plot(x_l, y_l, 'Color', color, 'LineWidth', 1);
    axis equal
    yticks(-0.12:0.04:0.12);
end

function x_d = discretizer(x, cpoint, gap_smallx, gap_largex)
    x_d = [];
    for i=1:length(x)
        if x(i)<=cpoint && ~mod(i-1, gap_smallx)
            x_d = [x_d, x(i)];
        elseif x(i)>cpoint && ~mod(i, gap_largex)
            x_d = [x_d, x(i)];
        end
    end
end

%%%%%%%%%%%%%%%%%%y_c y_t and theta calculator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y_c = get_y_c(x, m, k1)
    if x>=0 && x<=m
        y_c = k1/6*(x^3 - 3*m*x^2 + m^2*(3-m)*x);
    elseif x>m && x<=1
        y_c = k1/6*m^3*(1-x);
    else
        y_c = 0;
    end
end

function y_t = get_y_t(x, t)
    y_t = t/0.2*(0.2969*sqrt(x) - 0.126*x - 0.3516*x^2 +...
          0.2843*x^3 - 0.1036*x^4);
end

function theta = get_theta(x, m, k1)
    if x>=0 && x<=m
        theta = atan( k1/6*(3*x^2 - 6*m*x + m^2*(3-m)) );
    elseif x>m && x<=1
        theta = atan( -k1*m^3/6 );
    else
        theta = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%x_u/l and y_u/l calculator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function x_u = get_x_u(x, m, k1, t)
    x_u = x -  get_y_t(x, t)*sin( get_theta(x, m, k1));
end

function x_l = get_x_l(x, m, k1, t)
    x_l = x +  get_y_t(x, t)*sin( get_theta(x, m, k1));
end

function y_u = get_y_u(x, m, k1, t)
    y_u =  get_y_c(x, m, k1) +  get_y_t(x, t)*cos( get_theta(x, m, k1));
end

function y_l = get_y_l(x, m, k1, t)
    y_l =  get_y_c(x, m, k1) -  get_y_t(x, t)*cos( get_theta(x, m, k1));
end

function a_c = centerize(a)
    a_c = ([a(2:length(a)), 0] + a)/2;
    a_c = a_c(1:length(a_c)-1);
end