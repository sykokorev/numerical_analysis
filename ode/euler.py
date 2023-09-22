def euler(func, t_span, y0, h, args=()) -> float:

    '''
        Solves ODE IVP by using (explicit) Euler method
        -----------------------------------------------
        Paramters:      func: callable
                            Right-hand side of the equation. 
                            The time derivative of the state y at time t
                            func(t, x, *args)
                        t_span: 2-members sequence
                            Interval of integration (t0, tf)
                        y0: float
                            Initial state
                        steps: int
                            Number of integration steps
                        args: tuple, optional
                            Extra arguments to pass to the function
    '''

    ysol = [y0]
    tsol = [t_span[0]]
    n = (int)((t_span[1] - t_span[0]) / h)
    
    for i in range(n + 1):
        if tsol[i] + h <= t_span[1]:
            ysol.append(ysol[i] + h * func(tsol[i], ysol[i], *args))
            tsol.append(tsol[i] + h)

    if tsol[-1] != t_span[1]:
        h = t_span[1] - tsol[-1]
        ysol.append(ysol[-1] + h * func(tsol[-1], ysol[-1], *args))
        tsol.append(t_span[1])

    return tsol, ysol
