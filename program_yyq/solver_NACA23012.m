%NACA23012 vortex solver
global unit tangent normal deltas number v_infty alpha

x_l_c = flip( readmatrix('coord/x_l_c.csv') );
y_l_c = flip( readmatrix('coord/y_l_c.csv') );
x_l = flip( readmatrix('coord/x_l.csv') );
y_l = flip( readmatrix('coord/y_l.csv') );
x_u_c = readmatrix('coord/x_u_c.csv');
y_u_c = readmatrix('coord/y_u_c.csv');
x_u = readmatrix('coord/x_u.csv');
y_u = readmatrix('coord/y_u.csv');

unit = [[x_l_c.',y_l_c.'];[x_u_c.',y_u_c.']];
number = length(unit);
tangent_u = get_tangent(x_u,y_u);
tangent_l = get_tangent(x_l,y_l);
tangent = [tangent_l; tangent_u]; %'tangent' hasn't been and need not to be normalized
deltas = sqrt(tangent(:,1).^2 + tangent(:,2).^2);
normal = get_normal(tangent); %while 'normal' has been normalized
alpha = 5;
v_infty = 100;

main();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%main funciton%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function main()
    global unit tangent normal deltas number v_infty alpha
    gamma =  solver(number, unit, tangent, normal, v_infty, deg2rad(alpha), deltas);
    cl_zks = get_cl_zks(gamma, deltas, v_infty);
    writematrix(gamma, 'coord/gamma.csv');
    writematrix(unit, 'coord/unit.csv');
    writematrix(tangent, 'coord/tangent.csv');
    writematrix(deltas, 'coord/deltas.csv');
    writematrix(normal, 'coord/normal.csv');
    fl = fopen('coord/cl.txt', 'w');
    fprintf(fl, 'alpha = %d degree, cl_by_zhukovsky = %.4f\n', alpha, cl_zks);
    fclose(fl);
end

%%%%%%%%%%%%%%%%%%%%%%%%%get tangent and normal%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get_tangent
function tangent = get_tangent(x_p, y_p)
    point = [x_p.',y_p.'];
    tangent = ([point(2:length(point),:); [0,0]] - point);
    tangent = tangent(1:length(tangent)-1,:);
end

%get normal
function normal = get_normal(tangent)
    normal = zeros(length(tangent),2);
    for i=1:length(tangent)
        normal(i,:) = [-tangent(i,2), tangent(i,1)];
        normal(i,:) = normal(i,:)/sqrt(sum(normal(i,:).^2));
    end
end

%%%%%%%%%%%%%%%%%get velocity function and solve gamma%%%%%%%%%%%%%%%%%%%%%

%get velocity derived by stream function
function v = velocity(target, origin)
    d_square = sum((target-origin).^2);
    v_x = (target(2)-origin(2))/d_square;
    v_y = -(target(1)-origin(1))/d_square;
    v = [v_x, v_y];
end

%get matrix A for v_n, then we have A*v_n = b
function A_n = get_An(number, unit, tangent, normal, deltas)
    A_n = zeros(number);
    for i=1:number
        for j=1:number
            if j ~= i
                f = @(xi) dot( velocity(unit(i,:), unit(j,:)+tangent(j,:)*xi), normal(i,:) )*deltas(j);
                A_n(i,j) = integral(f, -1/2, 1/2, 'ArrayValued',true)/(2*pi);
            end
        end
    end
end

%solve the lambda
function gamma = solver(number, unit, tangent, normal, v_infty, alpha, deltas)
    A = get_An(number, unit, tangent, normal, deltas);
    b = zeros(number,1);
    vb = v_infty*[cos(alpha), sin(alpha)];

    A(number,:) = zeros(1,number);
    A(number,1) = 1;
    A(number,number) = 1; %kutta condition
    for i=1:number
        b(i) = -dot(vb, normal(i,:));
    end
    b(number) = 0;

    [P,R,C] = equilibrate(A);
    B = R*P*A*C;
    d = R*P*b;
    gamma = C*(B\d);
end

%%%%%%%%%%%%%%%%%%%%%%get cl by zhukovsky's method%%%%%%%%%%%%%%%%%%%%%%%%%

%get c_l
function c_l = get_cl_zks(gamma, deltas, v_infty)
    vortex = sum(gamma.*deltas);
    c_l = 2*vortex/v_infty;
end