using homeostasis

using Cyton

using DataFrames, Gadfly
# Gadfly defaults
Gadfly.set_default_plot_size(20cm, 20cm)
Gadfly.push_theme(Theme(background_color="white"))

function run(model::CellPopulation, runDuration::Time)
  print("Time to run:")
  @time begin
    counts = DataFrame(time=Time[], count=Int[])
    Δt = modelTimeStep(model)

    for tm in 1:Δt:runDuration
      Cyton.step(model)
      push!(counts, (tm, length(model.cells)))
    end
  end

  h = plot(counts, x=:time, y=:count, Geom.line)
  display(h)

  println("Done at model time=$(modelTime(model))")
end


@info rpad(lpad(" start ", 30, "-"), 55, "-")
@info "Creating population"
environment=environmentFactory()
population = createPopulation(100, (birth) -> homeostasis.cellFactory(birth),environment)
@info "Running model"
run(population, 100.0)
@info "done"
@info rpad(lpad(" end ", 31, "-"), 56, "-")

