import Cairo, Fontconfig
#
using Cyton

push!(LOAD_PATH, pwd())
using homeostasis

using DataFrames, Gadfly, Statistics, Serialization, ArgParse

# Gadfly defaults
Gadfly.set_default_plot_size(20cm, 20cm)
Gadfly.push_theme(Theme(
  background_color = "white", 
  alphas = [0.5],
  minor_label_font_size = 12pt,
  major_label_font_size = 14pt,
  point_label_font_size = 12pt,
  key_label_font_size = 12pt,
  key_title_font_size = 14pt,
  line_width = 2pt
  ))


  OUTPUT_DIR = "outputs"

#----------------------- Results -----------------------
struct RunResults
  parms::TrialParameters
  model::CytonModel
  generationCounts::DataFrame
  cellPhases::DataFrame
  proteinLevels::DataFrame
  divisions::Vector{Time}
end

function run(model::CytonModel, runDuration::Time, parms::TrialParameters)
  generationCounts = DataFrame(time=Time[], total=[], gen0=[], gen1=[], gen2=[], gen3=[], gen4=[], gen5=[], gen6=[], gen7=[], gen8=[], genOther=[])
  cellPhases = DataFrame(time=Time[], total=[], G1=[], S=[], G2M=[])
  proteinLevels = DataFrame(time=Time[], avgMyc=[], avgSurvivalProtein=[], il7=Float64[])

  # Record divisions in the model by registering an observer for Division events
  divisions = Time[]
  sizehint!(divisions, 100000)
  function divisionCounter(e::CellEvent, time::Time)
    if e == Division()
      push!(divisions, time)
    end
  end
  push!(model.eventCallbacks, divisionCounter)

  Δt = modelTimeStep(model)
  local tm
  for tm in 0:Δt:runDuration

    Cyton.step(model)
    cells = keys(model.cells)
    nCells = length(cells)

    if tm % 100 == 0
      @show tm
      @show nCells
    end

    if nCells == 0
      @warn "Stopping early - no life"
      break
    end

    local phaseCnts = Dict(G1 => 0, S => 0, G2M => 0)
    for cell in cells
      (p, _) = phase(cell, tm)
      phaseCnts[p] += 1
    end
    push!(cellPhases, (tm, length(cells), phaseCnts[G1], phaseCnts[S], phaseCnts[G2M]))

    local genCnts = zeros(10)
    for cell in cells
      gen = cell.generation
      if gen <=  8
        genCnts[gen+1] +=  1
      else
        genCnts[10] +=  1
      end
    end
    push!(generationCounts, (tm, length(cells), genCnts...))

    currentMyc = zeros(nCells)
    currentSp = zeros(nCells)
    for (i, cell) in enumerate(cells)
      myctimer = divisionTimer(cell)
      currentMyc[i] = myctimer.myc
      survivaltimer = survivalTimer(cell)
      currentSp[i] = survivaltimer.survivalProtein
    end
    if length(cells) == 0
      avgMyc = 0.0
      avgSp = 0.0
    else
      avgMyc = mean(currentMyc)
      avgSp = mean(currentSp)
    end
    push!(proteinLevels, (tm, avgMyc, avgSp, model.environmentAgents[1].concentration))
  end

  println("Done at model time = $(modelTime(model))")
  return RunResults(parms, model, generationCounts, cellPhases, proteinLevels, divisions)
end

function doRun(pSet)

  resultFileName = OUTPUT_DIR * "/results $(pSet).dat"

  if isfile(resultFileName)
    @info "$pSet already done, skipping run"
  else
    @info "-"^30 * " start trial $pSet " * "-"^30
    
    @info "Creating population"
    environment = environmentFactory()
    population = createPopulation(pSet.initPopulationSize, (birth) -> cellFactory(birth, pSet), environment)
    
    @info "Running model"
    result = run(population, 2000.0, pSet)
    Base.Filesystem.mkpath(OUTPUT_DIR)
    open(resultFileName, "w") do io
      serialize(io, result)
    end

  end

  plotResult(result)

  @info "-"^30 * " end trial $pSet " * "-"^30

end

function loadResults(parameters)
  results = RunResults[]
  for pSet in parameters
    fn = OUTPUT_DIR * "/results $(pSet).dat"
    if isfile(fn)
      result = deserialize(fn)
      push!(results, result)
    else
      @warn "$fn not found"
    end
  end

  return results
end

function plotResult(result)
  Base.Filesystem.mkpath(OUTPUT_DIR)

  vars = [:total :gen0 :gen1 :gen2 :gen3 :gen4 :gen5 :gen6 :gen7 :gen8 :genOther]
  title = "Cell counts=> $(result.parms)"
  h = plot(result.generationCounts, x=:time, y=Col.value(vars...), color=Col.index(vars...), Geom.line, Guide.title(title))
  fn = OUTPUT_DIR * "/" * title * ".png"
  h |> PNG(fn, 15cm, 15cm)

  vars = [:total :G1 :S :G2M]
  title = "Cell phases=> $(result.parms)"
  h = plot(result.cellPhases, x=:time, y=Col.value(vars...), color=Col.index(vars...), Geom.line, Guide.title(title))
  fn = OUTPUT_DIR * "/" * title * ".png"
  h |> PNG(fn, 15cm, 15cm)

  vars = [:avgMyc, :avgSurvivalProtein, :il7]
  title = "Protein=> $(result.parms)"
  h = plot(result.proteinLevels, x=:time, y=Col.value(vars...), color=Col.index(vars...), Geom.line, Guide.title(title))
  fn = OUTPUT_DIR * "/" * title * ".png"
  h |> PNG(fn, 15cm, 15cm)

  title = "Divisions=> $(result.parms)"
  h = plot(x=result.divisions, Geom.histogram, Guide.title(title), Guide.xlabel("time"))
  fn = OUTPUT_DIR * "/" * title * ".png"
  h |> PNG(fn, 15cm, 15cm)
end

function doPlots(results)
  Threads.@threads for result in results
    plotResult(result)
  end
end

function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table! s begin
      "--trial"
      help = "A trial number (to distinguish it from repeats with otherwise equal parameters)"
      arg_type = Int
      required = true

      "--initial-population"
      help = "Initial population size"
      arg_type = Int
      required = true

      "--threshold"
      help = "Threshold for the survival protein"
      arg_type = Float64
      required = true
    end

  return parse_args(s)
end

function main()
  parsed_args = parse_commandline()
  initPopSize = parsed_args["initial-population"]
  threshold = parsed_args["threshold"]
  trial = parsed_args["trial"]
  @show initPopSize
  @show threshold
  @show trial

  parms = TrialParameters(trial, threshold, initPopSize)
  doRun(parms)
end


main()
