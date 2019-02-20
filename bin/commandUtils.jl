using Pkg

function importSim(simFile)
  simName = split(basename(simFile), ".")[1]
  if (!startswith(simFile, "/"))
    simFile = "$(pwd())/$simFile"
  end
  include(simFile)
  simName
end

function activatePkg()
  Pkg.activate(abspath(PROGRAM_FILE, "..", ".."))
end
