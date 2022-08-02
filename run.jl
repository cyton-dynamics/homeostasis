import Cairo,Fontconfig
#
using Cyton
using homeostasis

using DataFrames, Gadfly, Plots, Statistics
# Gadfly defaults
Gadfly.set_default_plot_size(20cm, 20cm)
Gadfly.push_theme(Theme(
  background_color="white", 
  alphas=[0.5],
  minor_label_font_size = 12pt,
  major_label_font_size = 14pt,
  point_label_font_size = 12pt,
  key_label_font_size = 12pt,
  key_title_font_size = 14pt,
  line_width = 2pt
  ))

function run(model::CytonModel, runDuration::Time)
  print("Time to run:")
  @time begin
    counts = DataFrame(time=Time[], count=Int[])
    il7=DataFrame(time=Time[], il7=Float64[])
    generation = DataFrame(time=Time[], 
    total=[], 
    gen0 = [],
    gen1 = [],
    gen2 = [],
    gen3 = [],
    gen4 = [],
    gen5 = [],
    gen6 = [],
    gen7 = [],
    gen8 = [],
    genOther = []
    )

    cellphases = DataFrame(timer=Time[], 
    total=[], 
    G1 = [],
    S = [],
    G2M = []
    )

    myclevels=DataFrame(timer=Time[],
    avgMyc=[])


    timeInG1 = Dict{Cell,Int64}()
    


    Δt = modelTimeStep(model)

    for tm in 1:Δt:runDuration

      Cyton.step(model)
      push!(counts, (tm, length(model.cells)))
      push!(il7,(tm,model.environmentAgents[1].concentration))

      local phaseCnts = zeros(3)
      cells = model.cells
      for (cell,id) in cells
        divtimer = filter(x->typeof(x)==DivisionTimer,cell.timers)
        (phase,δt)=homeostasis.phase(divtimer[1])
        if phase <= 2
          phaseCnts[phase+1] += 1
        end
      end
      push!(cellphases, (tm, length(cells), phaseCnts...))

      local genCnts = zeros(10)
      cells = model.cells
      for (cell,id) in cells
        gen = cell.generation
        if gen <= 8
          genCnts[gen+1] += 1
        else
          genCnts[10] += 1
        end
      end
      push!(generation, (tm, length(cells), genCnts...))

      cells=model.cells
      currentMyc=[]
      for (cell,id) in cells
        myctimer=filter(x->typeof(x)==MycTimer,cell.timers)
        myctimer=myctimer[1]
        cellMyc=myctimer.myc
        push!(currentMyc,cellMyc)
      end
      avgMyc=mean(currentMyc)
      push!(myclevels,(tm,avgMyc))
        
    end
    cells = model.cells
    for (cell,id) in cells
      divtimer = filter(x->typeof(x)==DivisionTimer,cell.timers)
      divtimer=divtimer[1]
      timeSpentG1=divtimer.timeInG1
      push!(timeInG1,cell=>timeSpentG1)
    end
    dead=model.deadCells
    for (cell,id) in dead
      divtimer = filter(x->typeof(x)==DivisionTimer,cell.timers)
      divtimer=divtimer[1]
      timeSpentG1=divtimer.timeInG1
      push!(timeInG1,cell=>timeSpentG1)
    end


  end



  h = Gadfly.plot(counts, x=:time, y=:count, Geom.line, Coord.cartesian(xmin=0, xmax=10000, ymin=0, ymax=1500))
  #draw(PDF("CellCounts.pdf"), h)

  # h1 = plot(il7, x=:time, y=:il7, Geom.line)
  #draw(PDF("IL7.pdf"), h1)

  # vars2 = [:total :G1 :S :G2M ]
  # h3 = plot(cellphases, x=:timer, y=Col.value(vars2...), color=Col.index(vars2...),Geom.line)
  #draw(SVG("cellphases.svg"), h3)

  # vars = [:total :gen0 :gen1 :gen2 :gen3 :gen4 :gen5 :gen6 :gen7 :gen8 :genOther]
  # h2 = plot(generation, x=:time, y=Col.value(vars...), color=Col.index(vars...),Geom.line)
  #draw(PDF("generations.pdf"), h2)
  

  println("Done at model time=$(modelTime(model))")
return (h,cellphases,timeInG1,myclevels)
end


# @info rpad(lpad(" start ", 30, "-"), 55, "-")
# @info "Creating population"
# environment=environmentFactory()
# population = createPopulation(1000, (birth) -> cellFactory(birth),environment)
#population = createPopulation(50, (birth) -> cellFactory(birth))
# @info "Running model"
# run(population, 6000.0)
# @info "done"
# @info rpad(lpad(" end ", 31, "-"), 56, "-")



struct RunResults
  plot::Union{Nothing, Plot}
  phaseCounts
  timeInG1:: Dict{Cell,Int64}
  myclevels
end



parameters=trialParameters[]
results=RunResults[]
thresholdList=[log(1.1),1.1*log(1.1),1.2*log(1.1),1.3*log(1.1),1.4*log(1.1),1.5*log(1.1),1.6*log(1.1)]


for threshold in thresholdList
  currentSet=trialParameters(threshold)
  push!(parameters,currentSet)
end

for (index,pSet) in enumerate(parameters)
  local model,rt,phaseCounts,timeInG1,result
  @info rpad(lpad(" start trial $index ", 30, "-"), 55, "-")
  @info "Creating population"
  environment=environmentFactory()
  population = createPopulation(1000, (birth) -> cellFactory(birth,pSet),environment)
  @info "Running model"
  (plot,phaseCounts,timeInG1,myclevels)=run(population,3500.0)
  result=RunResults(plot,phaseCounts,timeInG1,myclevels)
  push!(results,result)
  @info "done"
  @info rpad(lpad(" end trial $index", 31, "-"), 56, "-")
end


@info "Running done. Summarising the data"

countSummary=DataFrame(
  Thresholds=[],
  counts=[],
  G1 = [],
  S = [],
  G2M = []
)

global MycSummary=DataFrame(
  timer=[],
  avgMyc=[]
)
global allThresholdMyc=[]

global allThresholdG1=[]
global timesInG1=DataFrame( timeInG1=[])



for (index,result) in enumerate(results)

  local cellphases=result.phaseCounts
  cellphases=Matrix(cellphases)
  local rowlength=size(cellphases,1)
  local phaseCnts=cellphases[rowlength-100:rowlength,:]
  phaseCnts=vec(mean(phaseCnts, dims=1))
  popfirst!(phaseCnts)
  push!(countSummary,(thresholdList[index],phaseCnts...))

  
  global allThresholdMyc=vcat(allThresholdMyc,thresholdList[index]*ones(size(result.myclevels,1)))
  global MycSummary=vcat(MycSummary,result.myclevels)


  local age=DataFrame(timeInG1=collect(values(result.timeInG1)))
  global allThresholdG1=vcat(allThresholdG1,thresholdList[index]*ones(size(age,1)))
  global timesInG1=vcat(timesInG1,age)




end
mycDF=hcat(allThresholdMyc,MycSummary)
timeInG1DF=hcat(allThresholdG1,timesInG1)

@info "Done"






# @info "Constructing cellphases plot"
# p1=groupedbar([countSummary[:,:G1] countSummary[:,:S] countSummary[:,:G2M]],
#   bar_position = :stack,
#   bar_width=0.4,
#   xticks=(1:length(thresholdList),round.(countSummary[:,:Thresholds],sigdigits=4)),
#   label=["G1" "S" "G2M"])

# p1=Plots.plot(countSummary[:,:Thresholds],[countSummary[:,:counts] countSummary[:,:G1] countSummary[:,:S] countSummary[:,:G2M]] ,
#     xlabel="Thresholds",
#     ylabel="Cell Counts",
#     label=["Total" "G1" "S" "G2M"],
#     leg=:bottomleft,
#     lw = 3,
#     marker = ([:hex :d :d :d])
# )

# timeInG1DF |>
# @vlplot(
#     :rect,
#     width=300, height=200,
#     x={:timeInG1, bin={maxbins=60}},
#     y={:x1, bin={maxbins=6},title="Thresholds"},
#     color="count()",
#     config={
#         range={
#             heatmap={
#                 scheme="greenblue"
#             }
#         },
#         view={
#             stroke="transparent"
#         }
#     }
# )|> save("timeinG1.pdf")


# mycDF |>
# @vlplot(
#     :line,
#     width=600, height=400,
#     x={:timer,axis={tickCount=3}},
#     y={:avgMyc,scale = {type="log"}},
#     color=:x1
# )|> save("AvgMyc.pdf")

