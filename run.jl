using Cyton
import Cyton: shouldDie, shouldDivide, inherit, step, stimulate, FateTimer, Stimulus


using homeostasis
import homeostasis: interact,environmentFactory,IL7,DeathTimer,shouldDivide, cellFactory,DivisionDestiny,MycTimer,step,inherit,destinyReached, update, DproTimer,DivisionTimer




using DataFrames, Gadfly
# Gadfly defaults
Gadfly.set_default_plot_size(20cm, 20cm)
Gadfly.push_theme(Theme(background_color="white"))

function run(model::CytonModel, runDuration::Time)
  print("Time to run:")
  @time begin
    counts = DataFrame(time=Time[], count=Int[])
    il7=DataFrame(time=Time[], il7=Float64[])
    Δt = modelTimeStep(model)

    for tm in 1:Δt:runDuration
      Cyton.step(model)
      push!(counts, (tm, length(model.cells)))
      push!(il7,(tm,model.environmentAgents[1].concentration))
    end
  end

  h = plot(counts, x=:time, y=:count, Geom.line)
  display(h)
  # h1 = plot(il7, x=:time, y=:il7, Geom.line)
  # display(h1)


  println("Done at model time=$(modelTime(model))")
end


@info rpad(lpad(" start ", 30, "-"), 55, "-")
@info "Creating population"
environment=environmentFactory()
population = createPopulation(25, (birth) -> cellFactory(birth),environment)
#population = createPopulation(50, (birth) -> cellFactory(birth))
@info "Running model"
run(population, 200.0)
@info "done"
@info rpad(lpad(" end ", 31, "-"), 56, "-")

