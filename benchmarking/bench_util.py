import os


def get_core_list():
    """

    """
    num_list = [1, 2, 4, 8, 16, 32, 64]
    num_cores = int(os.popen("sysctl -n hw.ncpu").read())  # macOS
    # num_cores = os.popen("nproc --all").read()  # Linux
    return sorted(i for i in num_list if i <= num_cores) 
