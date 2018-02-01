# Writes a makeflow file containing all necessary commands to run the 1D Sod
# convergence tests

# default configuration options
configure_options = {"disable-vec": "no option",
                     "with-hydro": "gizmo",
                     "with-riemann-solver": "exact",
                     "with-equation-of-state": "ideal-gas",
                     "with-hydro-dimension": 3,
                     "with-ext-potential": "none",
                     "enable-mpi": "no",
                     "enable-debug": "no",
                     "enable-debugging-checks": "no",
                     "enable-optimization": "yes",
                     "enable-sanitizer": "yes"}

# parameters: resolutions to run, and schemes to include
#resolutions = [100, 200, 400, 800, 1600, 3200]
resolutions = [100]
schemes = ["gizmo", "gadget2", "hopkins"]
# run command for a single simulation
run_command = "../swift -s -t 1 sodShock.yml 2>&1 | tee sodShock.log"

# get the compilation command to configure and compile the code for the given
# branch, number of dimensions and hydro scheme, using the given number of
# threads during compilation
def get_compilation_command(ndim, scheme, branch, nthread):
  global configure_options

  # git
  git_command = "git clone https://gitlab.cosma.dur.ac.uk/swift/swiftsim.git"
  # note that we move into the swiftsim folder here
  git_command += "; cd swiftsim"
  # we only switch branch if requested
  if not branch == "master":
    git_command += "; git checkout " + branch

  # configuration
  configure_options["with-hydro-dimension"] = ndim
  configure_options["with-hydro"] = scheme
  configure_command = "./autogen.sh"
  configure_command += "; ./configure"
  for configure_option in configure_options:
    configure_command += " --" + configure_option
    if not configure_options[configure_option] == "no option":
      configure_command += "={0}".format(configure_options[configure_option])

  # compilation
  # we move back to the original folder at the end
  make_command = "make -j {0}; cd ..".format(nthread)

  return git_command + "; " + configure_command + "; " + make_command

# create the makeflow file
wfile = open("SodShock_1D.makeflow", "w")
# resource information
wfile.write("CATEGORY=\"simulation\"\n")
wfile.write("CORES=1\n")
wfile.write("DISK=1000\n")
wfile.write("MEMORY=1000\n\n")

# commands to generate the initial conditions
# these are independent of the hydro scheme
for n in resolutions:
    wfile.write("ic_{0}.hdf5->sodShock.hdf5: makeIC.py\n".format(n))
    wfile.write("\tpython makeIC.py {0}\n\n".format(n))

# commands uses to run and analyze the simulations
for scheme in schemes:
  compilation = get_compilation_command(1, scheme, "master", 1)
  for n in resolutions:
    # we create a new folder to run the test
    folder_short = "swiftsim/examples/this_test"

    # input and output for the simulation
    input = "sodShock.yml->sodShock.yml".format(folder_short)
    input += " ic_{0}.hdf5->sodShock.hdf5".format(n, folder_short)
    output = "lastsnap_{0}_{1}.hdf5->{2}/sodShock_0001.hdf5".format(
      scheme, n, folder_short)
    output += " timesteps_{0}_{1}.txt->{2}/timesteps_1.txt".format(
      scheme, n, folder_short)
    output += " log_{0}_{1}.log->{2}/sodShock.log".format(
      scheme, n, folder_short)

    # prepare the simulation directory
    # we need to create the new folder and manually move all input files
    # after this is done, we make the new folder the working directory
    setup_folder = "mkdir {0}".format(folder_short)
    setup_folder += "; mv sodShock.yml {0}/".format(folder_short)
    setup_folder += "; mv sodShock.hdf5 {0}/".format(folder_short)
    setup_folder += "; cd {0}".format(folder_short)

    # now generate the simulation command
    wfile.write("{0}: {1}\n".format(output, input))
    wfile.write("\t{0}; {1}; {2}\n\n".format(
      compilation, setup_folder, run_command))

    # simulation analysis
    # this is done separately
    input = "lastsnap_{0}_{1}.hdf5->sodShock_0001.hdf5".format(scheme, n)
    input += " timesteps_{0}_{1}.txt->timesteps_1.txt".format(scheme, n)
    input += " analyze.py"
    output = "summary_{0}_{1}.txt->summary.txt".format(scheme, n)
    output += " result_{0}_{1}.png->result.png".format(scheme, n)
    wfile.write("{0}: {1}\n".format(output, input))
    wfile.write("\tpython analyze.py\n\n")
