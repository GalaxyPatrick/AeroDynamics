%for calculating cp of NACA23012

global unit tangent deltas normal number v_infty alpha gamma x_d_c

gamma = readmatrix('coord/gamma.csv');
unit = readmatrix('coord/unit.csv');
tangent = readmatrix('coord/tangent.csv');
deltas = readmatrix('coord/deltas.csv');
normal = readmatrix('coord/normal.csv');
x_d_c = readmatrix('coord/x_d_c.csv');
number = length(unit);
v_infty = 100;
alpha = 5;

main()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%main funciton%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main()
    global unit tangent deltas normal number v_infty alpha gamma x_d_c
    cp = get_cp(unit, gamma, tangent, v_infty, number, deg2rad(alpha), deltas);
    cl_intp = get_cl_intp(cp, deltas, number, normal);
    cp = cp.';
    scatter( x_d_c, flip(cp(1:number/2)) ,18 ,[1,0,0]);
    hold on
    scatter( x_d_c, cp(number/2+1:number) ,18 ,[0,0,1]);
    fl = fopen('coord/cl.txt', 'a');
    fprintf(fl, 'alpha = %d degree, cl_by_pressure_integral = %.4f\n', alpha, cl_intp);
    fclose(fl);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%calculate cp%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = velocity(target, origin)
    d_square = sum((target-origin).^2);
    v_x = (target(2)-origin(2))/d_square;
    v_y = -(target(1)-origin(1))/d_square;
    v = [v_x, v_y];
end

function [c_p, c_l] = get_cp(unit, gamma, tangent, v_infty, number, alpha, deltas)
    c_p = zeros(number,1);
    v_t = zeros(number,1);
    c_l = 0;
    tan_norm = zeros(number,2);
    for i=1:length(tangent)
        tan_norm(i,:) = tangent(i,:)/deltas(i);
    end

    for i=1:number
        mid = zeros(number,1);
        for j=1:number
            if j==i
                mid(j) = gamma(j)/2;
            else
                f = @(xi) dot( velocity(unit(i,:), unit(j,:)+tangent(j,:)*xi), tan_norm(i,:) )*deltas(j);
                mid(j) = gamma(j)*integral(f, -1/2, 1/2, 'ArrayValued',true)/(2*pi);
            end
        end
        v_t(i) = sum(mid) + dot(v_infty*[cos(alpha), sin(alpha)], tan_norm(i,:));
        c_p(i) = 1-(v_t(i)/v_infty)^2;
    end
end

function c_l = get_cl_intp(c_p, deltas, number, normal)
    c_l = 0;
    for i=1:number
        c_l = c_l - c_p(i)*deltas(i)*normal(i,2);
    end
end