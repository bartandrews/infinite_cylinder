import os
import sys


if __name__ == '__main__':

    for path in sys.path:
        print(path)

    print("sys.executable = ", sys.executable)

    os.system("nohup " + sys.executable + " code/ground_state.py -thr 1 -mod 'FerHofSqu1' -chi 500 -t1 1 -V 10 -Vtype 'Coulomb' -Vrange 1 -n 1 9 -nphi 1 3 -LxMUC 1 -Ly 9")
