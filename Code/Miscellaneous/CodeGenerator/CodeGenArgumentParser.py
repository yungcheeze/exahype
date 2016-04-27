import os
import AvailableConfigs

def validateLibxsmmGenerator(i_parser, i_arg):
    if not os.path.isdir(i_arg):
        i_parser.error("The libxsmm directory {} does not exist".format(i_arg))
    else:
        # directory exists, but is it a valid path to the libxsmm generator backend?
        l_pathToLibxsmm = i_arg
        l_pathToLibxsmmGenerator = l_pathToLibxsmm + "/bin/libxsmm_gemm_generator"
        if(os.path.isfile(l_pathToLibxsmmGenerator)):
            return l_pathToLibxsmm+"/bin"
        else:
            i_parser.error("Can't find the code generator of libxsmm. Did you run 'make generator' to build the library?")


def validateArchitecture(i_parser, i_arg):
    l_architecture = str(i_arg)

    # when we have to postprocess the generated assembly code we may
    # support only a subset of the available microarchitectures
    if(l_architecture not in AvailableConfigs.architectures):
        print("Driver: Unkown or unsupported microarchitecture. Continue with noarch")
        l_architecture = 'noarch'

    return l_architecture


def validatePrecision(i_parser, i_arg):
    l_precision = str(i_arg)

    if(l_precision not in AvailableConfigs.precisions):
        print("Unknown precision specified. Continue with double precision")
        l_precision    = 'DP' 

    return l_precision


def validateNumerics(i_parser, i_arg):
    l_numericsType = i_arg

    if(l_numericsType not in AvailableConfigs.numerics):
        i_parser.error("Numerics not supported. Available options are " + str(AvailableConfigs.numerics).format(i_arg))

    return l_numericsType
