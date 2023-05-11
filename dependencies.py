from os import environ as env
from pathlib import Path
from shutil import which

def get_executables(program):
    #checks if program's executable is reachable
    #and returns location or raises an error
    env_path_executable = check_executable_in_user_envs(program, verbose=False)
    if env_path_executable:
        executable = str(env_path_executable / program["executable"])
    elif check_is_executable_is_in_PATH(program["executable"], verbose=False):
        executable =  program["executable"]
    else:
        user_path = program["executable"]
        msg = "{} is not defined neither in"
        msg += " ${} nor $PATH"
        raise RuntimeError(msg.format(program["executable"], program["user_path"]))
    return executable


def check_executable_in_user_envs(user_env, program="default", verbose=True):
    #This function checks if user has defined a path
    #for program needed to run orgpba2
    executable = user_env["executable"]
    user_path = user_env["user_path"]
    #Try to find if variable is user_path is defined
    try:
        user_path_check = Path(env[user_path])
    except KeyError:
        if verbose:
            msg = "user env {} is not defined for {}"
            print(msg.format(user_path, program))
        return False
    #If user_path is defined, try to find all executables
    #needed to run the program
    else:
        full_path = Path(user_path_check) / executable
        check = full_path.absolute().exists()
        if not check:
            if verbose:
                msg = "{} doesn't contain {} executable"
                print(msg.format(user_path_check, executable))
            return False
        else:
            if verbose:
                print("{} found at {}".format(executable, full_path))
            return user_path_check


def check_is_executable_is_in_PATH(executable, verbose=True):
    #Checks if executable is included in $PATH variable
    msg = "{} is not defined in $PATH."
    msg += " {} path must be set either in $PATH"
    msg += " or set in user env path"
    check = which(executable)
    if check:
        if verbose:
            print("{} found at {}".format(executable, check))
    else:
        if verbose:
            print(msg.format(executable, executable))
        return False
    return True

