# Master file for PulseBlast2
# Sam Wallis-Riches, Kaine Bunting, 2019

from rfi_mitigation import SigmaClip_Mitigator, Bayesian_Mitigator, NN_Mitigator
import template_builder as tb
from flux_calibrator import FluxCalibrator

print("""
Welcome to PulseBlast V2
""")

# have verbose parsed into each 'run' function

errors = []

dirs = []
psr_names = []

verbose = False
r_verbose = False
c_verbose = False
m_verbose = False
t_verbose = False

frontend = None
subbands = None
template_dir = None

iterations = None
temp_dir = None
saveddata_dir = None
epoch_avg = None

cont_name = None
cont_dir = None
saveddata_dir2 = None

c_index = -1
r_index = -1
t_index = -1

commands = input("""Please type your command in here... """)
print("   ")
if commands == "help":
    print("""
Commands should be written in the format:

directories psr_names flags_and_parameters


directories                 Can be a single string, list or a .txt file
                            that contains a list of directories containing
                            .fits files.

psr_names                   Can be a single string, list or a .txt file
                            that contains a list of pulsar names.

functions_and_parameters    Any of the commands below.
                            They will run in the order they are given,
                            except template builder which will first if
                            its flag "-m" is given anywhere after directories
                            and psr_names


Running the code in verbose is an option:
Include the -v flag in your command within flags_and_parameters
eg. dirs.txt psr_names.txt -v -r

If you want only certain processes to run in verbose, you can set them to true
by adding the letters of the function tags after -v
eg. -v mc means run template builder and calibration in verbose

Verbose -v [m][c][r][t]
-v                          Verbose
m                           Template builder verbose
c                           Calibration verbose
r                           RFI mitigation verbose
t                           Timing verbose


Below, any optional commands are in square [] brackets.
The data type is given in curly {} brackets.
DO NOT type any parameters in "", just put the information after its flag,  unless in a list.
DO NOT use spaces if specifing directories/pulsar names in a list.
eg. ["directory1","diretory2","diretory3"]
If optional flags are not included, then the program will use defaults.

======================================================================================================
Template Builder: -m -f "frontened" -s subbands [-p "template directory"]
-m                          Template builder
-f "frontend"               Specifies frontend {str} (required)
-s subbands                 Specifiies subbands {int} (required)
-p "template directory"     Specificies template directory {str} (optional)

RFI Mitigation: -r/-rs/-rb/-rn [-i iterations] [-q "temp_dir"] [-g "saveddata_dir"] [-e]
-r                          RFI mitigation (defaults to mode S)
-rs                         RFI mitigation (mode S)
-rb                         RFI mitigation (mode B)
-rn                         RFI mitigation (mode N)
-i iterations               Specifies iterations {int} (optional)
-q "temp_dir"               Specifies template directory {str} (optional)
-g "saveddata_dir"          Specifies saved data directory {str} (optional)
-e                          Specifies epoch average to be True (defaults to false if not given)

Flux Calibration: -c -n "cont_name" [-w "cont_dir"] [-d "saveddata_dir"]
-c                          Flux calibration
-n "cont_name"              Specifies continuum name {str} (required)
-w "cont_dir"               Specifies continuum directory {str} (optional)
-d "saveddata_dir"          Specifies saved date directory {str} (optional)

Timing: -t
-t                          Calculate TOAs
======================================================================================================

""")

elif commands == "exit":
    exit()

else:
    commands = commands.split(" ")
    if ("[" in commands[0]) and ('"' in commands):
        dirs = commands[0]
    elif ".txt" in commands[0]:
        try:
            with open(commands[0], "r") as f:
                for line in f:
                    if line[-1] == "\n":
                        dirs.append(line[:-1])
                    else:
                        dirs.append(line)
        except FileNotFoundError:
            print("Directory .txt file not found")
            exit()
    else:
        dirs.append(commands[0])

    if "[" in commands[1]:
        psr_names = commands[1]
    elif ".txt" in commands[1]:
        try:
            with open(commands[1], "r") as f:
                for line in f:
                    if line[-1] == "\n":
                        psr_names.append(line[:-1])
                    else:
                        psr_names.append(line)
        except FileNotFoundError:
            print("Pulsar names .txt file not found")
            exit()
    else:
        psr_names.append(commands[1])

#testing
    print(commands)
    del commands[0]
    del commands[0]
    print(commands)

    if "-v" in commands:
        try:
            if "-" not in commands[commands.index("-v") + 1]:
                v_params = list(commands[commands.index("-v") + 1])
                if "m" in v_params:
                    m_verbose = True
                if "c" in v_params:
                    c_verbose = True
                if "r" in v_params:
                    r_verbose = True
                if "t" in v_params:
                    t_verbose = True
                del commands[commands.index("-v") + 1]
            elif "-" in commands[commands.index("-v") + 1]:
                r_verbose = True
                m_verbose = True
                c_verbose = True
                t_verbose = True
        except IndexError:
            r_verbose = True
            m_verbose = True
            c_verbose = True
            t_verbose = True
        del commands[commands.index("-v")]

    if "-m" in commands:
        print("first")
        if "-f" in commands:
            try:
                if "-" not in commands[commands.index("-f") + 1]:
                    try:
                        frontend = str(commands[commands.index("-f") + 1])
                    except ValueError:
                        error = "Invalid frontend given"
                        errors.append(error)
                elif "-" in commands[commands.index("-f") + 1]:
                    error = "Frontend not supplied"
                    errors.append(error)
            except IndexError:
                error = "Frontend not supplied"
                errors.append(error)
        elif "-f" not in commands:
            error = "Frontend not supplied"

        if "-s" in commands:
            try:
                if "-" not in commands[commands.index("-s") + 1]:
                    try:
                        subbands = int(commands[commands.index("-s") + 1])
                    except ValueError:
                        error = "Invalid subband number given"
                        errors.append(error)
                    if subbands <= 0:
                        subbands = None
                        error = "Invalid subband number given"
                        errors.append(error)
                elif "-" in commands[commands.index("-s") + 1]:
                    error = "Subband number not supplied"
                    errors.append(error)
            except IndexError:
                error = "Subband number not supplied"
                errors.append(error)
        elif "-s" not in commands:
            error = "Subband number not supplied"

        if "-p" in commands:
            try:
                if "-" not in commands[commands.index("-p") + 1]:
                    try:
                        template_dir = str(commands[commands.index("-p") + 1])
                    except ValueError:
                        error = "Invalid template directory given"
                        errors.append(error)
                elif "-" in commands[commands.index("-p") + 1]:
                    error = "Subband number not supplied"
                    errors.append(error)
            except IndexError:
                error = "Subband number not supplied"
                errors.append(error)

        if (frontend != None) and (subbands != None):
            if template_dir == None:
                # run without template_dir
                pass
            elif template_dir != None:
                # run with template_dir
                pass

    counter = 0
    for elem in commands:
        if "-r" in elem:
            counter += 1
    if counter == 1:
        if ("-r" in commands) or ("-rs" in commands) or ("-rn" in commands) or ("-rb" in commands):
            if "-i" in commands:
                try:
                    if "-" not in commands[commands.index("-i") + 1]:
                        try:
                            iterations = int(commands[commands.index("-i") + 1])
                        except ValueError:
                            error = "Invalid iterations number given"
                            errors.append(error)
                    elif "-" in commands[commands.index("-i") + 1]:
                        error = "Iterations not supplied"
                        errors.append(error)
                except IndexError:
                    error = "Iterations not supplied"
                    errors.append(error)

            if "-q" in commands:
                try:
                    if "-" not in commands[commands.index("-q") + 1]:
                        try:
                            temp_dir = str(commands[commands.index("-q") + 1])
                        except ValueError:
                            error = "Invalid template directory given"
                            errors.append(error)
                    elif "-" in commands[commands.index("-q") + 1]:
                        error = "Template diretory not supplied"
                        errors.append(error)
                except IndexError:
                    error = "Template diretory not supplied"
                    errors.append(error)

            if "-g" in commands:
                try:
                    if "-" not in commands[commands.index("-g") + 1]:
                        try:
                            saveddata_dir = str(commands[commands.index("-g") + 1])
                        except ValueError:
                            error = "Invalid save data diretory given"
                            errors.append(error)
                    elif "-" in commands[commands.index("-g") + 1]:
                        error = "Save data directory not supplied"
                        errors.append(error)
                except IndexError:
                    error = "Save data directory not supplied"
                    errors.append(error)

            if "-e" in commands:
                epoch_avg = True


    elif counter != 1 and counter != 0:
        error = "Too many RFI mode flags given"
        errors.append(error)

    if "-c" in commands:
        if "-n" in commands:
            try:
                if "-" not in commands[commands.index("-n") + 1]:
                    try:
                        cont_name = str(commands[commands.index("-n") + 1])
                    except ValueError:
                        error = "Invalid continuum name given"
                        errors.append(error)
                elif "-" in commands[commands.index("-n") + 1]:
                    error = "Continuum name not supplied"
                    errors.append(error)
            except IndexError:
                error = "Continuum name not supplied"
                errors.append(error)
        elif "-n" not in commands:
            error = "Continuum name not supplied"
            errors.append(error)

        if "-w" in commands:
            try:
                if "-" not in commands[commands.index("-w") + 1]:
                    try:
                        cont_dir = str(commands[commands.index("-w") + 1])
                    except ValueError:
                        error = "Continuum directory not given"
                        errors.append(error)
                elif "-" in commands[commands.index("-w") + 1]:
                    error = "Continuum diretory not supplied"
                    errors.append(error)
            except IndexError:
                error = "Continuum diretory not supplied"
                errors.append(error)
        elif "-w" not in commands:
            error = "Continuum diretory not supplied"
            errors.append(error)

        if "-d" in commands:
            try:
                if "-" not in commands[commands.index("-d") + 1]:
                    try:
                        saveddata_dir2 = str(commands[commands.index("-d") + 1])
                    except ValueError:
                        error = "Save data directory not given"
                        errors.append(error)
                elif "-" in commands[commands.index("-d") + 1]:
                    error = "Save data directory not supplied"
                    errors.append(error)
            except IndexError:
                error = "Save data directory not supplied"
                errors.append(error)

    if "-t" in commands:
        # import timing parameters from input
        pass

    if "-r" in commands:
        r_index = commands.index("-r")
    elif "-rs" in commands:
        r_index = commands.index("-rs")
    elif "-rb" in commands:
        r_index = commands.index("-rb")
    elif "-rn" in commands:
        r_index = commands.index("-rn")

    if "-c" in commands:
        c_index = commands.index("-c")

    if "-t" in commands:
        t_index = commands.index("-t")


# refer to previous operations' saveddata_dir not its own
    if c_index < r_index:
        # run Calibration first, then rfi
        if t_index < c_index:
            #run timing first
            pass

        if ("-c" in commands) and (cont_name != None) and (cont_dir != None):
            if saveddata_dir2 != None:
                #run cal with saveddata_dir2
                pass
            elif saveddata_dir2 == None:
                # run cal
                saveddata_dir = saveddata_dir2
                pass

        if (t_index > c_index) and (t_index < r_index):
            # run timing
            pass

        if (("-r" in commands) or ("-rs" in commands) or ("-rn" in commands) or ("-rb" in commands)) and (counter == 1):
            if epoch_avg == None:
                epoch_avg = False
            if iterations != None:
                if temp_dir != None:
                    if saveddata_dir != None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with iterations, temp_dir, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, temp_dir = temp_dir, saveddata_dir = saveddata_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with iterations, temp_dir, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, temp_dir = temp_dir, saveddata_dir = saveddata_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with iterations, temp_dir, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, temp_dir = temp_dir, saveddata_dir = saveddata_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                    elif saveddata_dir == None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with iterations, temp_dir, epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, temp_dir = temp_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with iterations, temp_dir, epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, temp_dir = temp_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with iterations, temp_dir, epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, temp_dir = temp_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )

                elif temp_dir == None:
                    if saveddata_dir != None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with iterations saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, saveddata_dir = saveddata_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with iterations, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, saveddata_dir = saveddata_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with iterations, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, saveddata_dir = saveddata_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                    elif saveddata_dir == None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with iterations, epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with iterations, epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with iterations, epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
            elif iterations == None:
                if temp_dir != None:
                    if saveddata_dir != None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with temp_dir, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, temp_dir = temp_dir, saveddata_dir = saveddata_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with temp_dir, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, temp_dir = temp_dir, saveddata_dir = saveddata_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with temp_dir, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, temp_dir = temp_dir, saveddata_dir = saveddata_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                    elif saveddata_dir == None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with temp_dir, epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, temp_dir = temp_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with temp_dir, epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, temp_dir = temp_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with temp_dir, epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, temp_dir = temp_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                elif temp_dir == None:
                    if saveddata_dir != None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, saveddata_dir = saveddata_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, saveddata_dir = saveddata_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, saveddata_dir = saveddata_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                    elif saveddata_dir == None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, epoch_avg = epoch_avg, verbose = r_verbose )

            rfi.mitigation_setup()

        if t_index > r_index:
            # run timing last
            pass

    elif r_index < c_index:
        # run rfi first, then cal
        if t_index < r_index:
            # run timing first
            pass

        if (("-r" in commands) or ("-rs" in commands) or ("-rn" in commands) or ("-rb" in commands)) and (counter == 1):
            if epoch_avg == None:
                epoch_avg = False
            if iterations != None:
                if temp_dir != None:
                    if saveddata_dir != None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with iterations, temp_dir, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, temp_dir = temp_dir, saveddata_dir = saveddata_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with iterations, temp_dir, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, temp_dir = temp_dir, saveddata_dir = saveddata_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with iterations, temp_dir, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, temp_dir = temp_dir, saveddata_dir = saveddata_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                    elif saveddata_dir == None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with iterations, temp_dir, epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, temp_dir = temp_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with iterations, temp_dir, epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, temp_dir = temp_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with iterations, temp_dir, epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, temp_dir = temp_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )

                elif temp_dir == None:
                    if saveddata_dir != None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with iterations saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, saveddata_dir = saveddata_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with iterations, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, saveddata_dir = saveddata_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with iterations, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, saveddata_dir = saveddata_dir, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                    elif saveddata_dir == None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with iterations, epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with iterations, epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with iterations, epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, iterations = iterations, epoch_avg = epoch_avg, verbose = r_verbose )
            elif iterations == None:
                if temp_dir != None:
                    if saveddata_dir != None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with temp_dir, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, temp_dir = temp_dir, saveddata_dir = saveddata_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with temp_dir, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, temp_dir = temp_dir, saveddata_dir = saveddata_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with temp_dir, saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, temp_dir = temp_dir, saveddata_dir = saveddata_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                    elif saveddata_dir == None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with temp_dir, epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, temp_dir = temp_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with temp_dir, epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, temp_dir = temp_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with temp_dir, epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, temp_dir = temp_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                elif temp_dir == None:
                    if saveddata_dir != None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, saveddata_dir = saveddata_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, saveddata_dir = saveddata_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with saveddata_dir and epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, saveddata_dir = saveddata_dir, epoch_avg = epoch_avg, verbose = r_verbose )
                    elif saveddata_dir == None:
                        if ("-r" in commands) or ("-rs" in commands):
                            # run sigma clip with epoch_avg=true
                            for name in psr_names:
                                rfi = SigmaClip_Mitigator( name, *dirs, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rb" in commands:
                            # run bayesian with epoch_avg=true
                            for name in psr_names:
                                rfi = Bayesian_Mitigator( name, *dirs, epoch_avg = epoch_avg, verbose = r_verbose )
                        elif "-rn" in commands:
                            # run neural with epoch_avg=true
                            for name in psr_names:
                                rfi = NN_Mitigator( name, *dirs, epoch_avg = epoch_avg, verbose = r_verbose )

            rfi.mitigation_setup()

        if (t_index > r_index) and (t_index < c_index):
            # run timing
            pass

        if ("-c" in commands) and (cont_name != None) and (cont_dir != None):
            print("Kaine")
            if saveddata_dir2 != None:
                # run cal with saveddata_dir2
                pass
            elif saveddata_dir2 == None:
                # run cal
                saveddata_dir2 = saveddata_dir
                pass

        if t_index > c_index:
            # run timing last
            pass

if ("-l" in commands) and (errors != []):
    # write to log file
    with open("errors.txt", "a") as f:
        for error in errors:
            f.write(error)
elif ("-l" not in commands) and (errors != []):
    print("""
If any errors occurred, the program continued with any
operations it could perform until it ran out.

A list of errors can be found below:
""")

    for error in errors:
        print(error)


# TESTING AREA ONLY:
print(dirs)
print(psr_names)
print(verbose)
print(r_verbose)
print(c_verbose)
print(m_verbose)
print(t_verbose)
print(frontend)
print(subbands)
print(template_dir)
print(iterations)
print(temp_dir)
print(saveddata_dir)
print(epoch_avg)
print(cont_name)
print(cont_dir)
print(saveddata_dir2)
print(c_index)
print(r_index)
