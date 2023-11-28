%yyq26314@outlook.com
%run all programs

global need_discretize need_solve_gamma need_calculate_cp
need_discretize = 0;
need_solve_gamma = 1;
need_calculate_cp = 1;

main()

function main()
    global need_discretize need_solve_gamma need_calculate_cp
    if need_discretize
        NACA23012
    end

    if need_solve_gamma
        solver_NACA23012
    end

    if need_calculate_cp
        cp_calculator
    end
end