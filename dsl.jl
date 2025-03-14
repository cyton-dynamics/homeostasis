
# using ModelingToolkit

# macro cell(ex)
#   println(cell)
# end

function parse(::LineNumberNode) end
function parse(stuff::Expr)
  for ex in stuff.args
    if ex isa LineNumberNode
      continue
    end

    head = ex.head
    args = ex.args
    if head != :(=)
      show(head)
    elseif  args[1] == :type
      println("type = $(ex.args[2])")
    end

    # println("$(ex.head) => $(ex.args)")
  end
end

macro shared(stuff)
  # parse(stuff)
end

macro cell(stuff)
  parse(stuff)
end

#D = Differential(t)

@shared begin
    type = "IL7"
    protein = IL7
    init = IL7 = 0
    time = D(IL7) = k
end


@cell begin
  type = "BaxBax KO"

  proteins  = [myc, survival, cell_cycle]
  events    = [death, division, cycle_start]
  states    = [G1 S G2M]
  callbacks = [cb1, cb2]

  on_init = begin 
    Tdeath = rand(λ_death_threshold)
    Tdiv   = rand(λ_division_threshold)
    state = out_of_cycle
  end

  update = begin
    D(myc) ~ o(IL7) - k*myc
    D(IL7) ~ -o(IL7)
    survival ~ f(t)
  end

  on_divide = begin
    myc -> myc + rand()
    survival -> survival
    state -> out_of_cycle
  end

  death = begin
    type = "death"
    trigger = survival < T
    action = die()
  end

  divide = begin
    type = "divide"
    trigger = state == in_cycle && myc > T 
    action = divide()
  end

  start_cycle = begin
    type = "start cycle"
    trigger = state != in_cycle && myc > T 
    action = state <- in_cycle
  end

  events = [death, divide, start_cycle]
end

