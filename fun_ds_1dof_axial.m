function dydt = fun_ds_1dof_axial(t,y)
    global m c k So
    if t <= 0.2
        Ff = So(1);
    elseif t > 0.2 && t <= 0.4
        Ff = So(2);
    elseif t > 0.4 && t <= 0.6
        Ff = So(3);
    elseif t > 0.6 && t <= 0.8
        Ff = So(4);
    else
        Ff = So(5);
    end
    dydt = [y(2); (1e1*Ff)/m - c/m*y(2) - k/m*y(1)];
