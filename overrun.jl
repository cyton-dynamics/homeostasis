using DifferentialEquations

test(u, t, integrator) = t-0.5
affect!(i) = add_tstop!(i, 2)

prob = ODEProblem((du, u, p, t) -> du[1]=0, [0], (0.0, 1.0))
sol = solve(prob, callback = ContinuousCallback(test, affect!))
sol.t[end]
