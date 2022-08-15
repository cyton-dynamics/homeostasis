import Cairo, Fontconfig
#
using Cyton

push!(LOAD_PATH, pwd())
using homeostasis

using DataFrames, Gadfly, Statistics, Serialization

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


#----------------------- Results -----------------------
struct RunResults
  parms::TrialParameters
  model::CytonModel
  generationCounts::DataFrame
  cellPhases::DataFrame
  proteinLevels::DataFrame
end

function run(model::CytonModel, runDuration::Time, parms::TrialParameters)
  generationCounts = DataFrame(time=Time[], total=[], gen0=[], gen1=[], gen2=[], gen3=[], gen4=[], gen5=[], gen6=[], gen7=[], gen8=[], genOther=[])
  cellPhases = DataFrame(time=Time[], total=[], G1=[], S=[], G2M=[])
  proteinLevels = DataFrame(time=Time[], avgMyc=[], avgSurvivalProtein=[], il7=Float64[])

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
  return RunResults(parms, model, generationCounts, cellPhases, proteinLevels)
end

function doRuns()
  parameters = TrialParameters[]
  thresholdList = [0.1, 0.2, 0.3]
  initPopSizes = [10, 100, 1000]
  nTrials = 5
  for threshold in thresholdList
    for initPopSize in initPopSizes
      for i in 1:nTrials
        currentSet = TrialParameters(i, threshold, initPopSize)
        push!(parameters, currentSet)
      end
    end
  end

  Threads.@threads for pSet in parameters
    @info "-"^30 * " start trial $pSet " * "-"^30
    
    @info "Creating population"
    environment = environmentFactory()
    population = createPopulation(pSet.initPopulationSize, (birth) -> cellFactory(birth, pSet), environment)
    
    @info "Running model"
    result = run(population, 2000.0, pSet)
    open("results $(pSet).dat", "w") do io
      serialize(io, result)
    end

    @info "-"^30 * " end trial $pSet " * "-"^30
  end

  return parameters
end

function loadResults(parameters)
  results = RunResults[]
  for pSet in parameters
    result = deserialize("results $(pSet).dat")
    push!(results, result)
  end

  return results
end

function doPlots(results)
  Base.Filesystem.mkpath("outputs")

  for result in results
    vars = [:total :gen0 :gen1 :gen2 :gen3 :gen4 :gen5 :gen6 :gen7 :gen8 :genOther]
    title = "Cell counts: $(result.parms)"
    h = plot(result.generationCounts, x=:time, y=Col.value(vars...), color=Col.index(vars...), Geom.line, Guide.title(title))
    fn = "outputs/" * title * ".png"
    h |> PNG(fn, 15cm, 15cm)

    vars = [:total :G1 :S :G2M]
    title = "Cell phases: $(result.parms)"
    h = plot(result.cellPhases, x=:time, y=Col.value(vars...), color=Col.index(vars...), Geom.line, Guide.title(title))
    fn = "outputs/" * title * ".png"
    h |> PNG(fn, 15cm, 15cm)

    vars = [:avgMyc, :avgSurvivalProtein, :il7]
    title = "Protein: $(result.parms)"
    h = plot(result.proteinLevels, x=:time, y=Col.value(vars...), color=Col.index(vars...), Geom.line, Guide.title(title))
    fn = "outputs/" * title * ".png"
    h |> PNG(fn, 15cm, 15cm)
  end
end

parameters = doRuns()
results = loadResults(parameters)
doPlots(results)