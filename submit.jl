using homeostasis
OUTPUT_DIR = "outputs"

function makeParameters()
  thresholdList = [0.1, 0.2, 0.3]
  initPopSizes = [10, 100, 1000]
  nTrials = 5
  
  parameters = TrialParameters[]
  for threshold in thresholdList
    for initPopSize in initPopSizes
      for i in 1:nTrials
        currentSet = TrialParameters(i, threshold, initPopSize)
        push!(parameters, currentSet)
      end
    end
  end

  return parameters  
end

parameters = makeParameters()

for parameter in parameters

  resultFileName = OUTPUT_DIR * "/results $(parameter).dat"
  if isfile(resultFileName)
    @info "$parameter already done, skipping run"
    continue
  end

  @info "Submitting $(parameter)"
  trial = parameter.trial
  threshold = parameter.threshold
  initPopSize = parameter.initPopulationSize
  cmd = `sbatch --cpus-per-task=2 --time=2-00:00:00 --mem=200GB --partition=regular --qos=bonus --job-name=il7-homeostasis batch.sh $(trial) $(threshold) $(initPopSize)`
  run(cmd)
end

@info "$(length(parameters)) jobs submitted"