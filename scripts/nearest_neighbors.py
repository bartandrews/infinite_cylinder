import numpy as np

if __name__ == '__main__':

    for i in range(1, 101):
        for j in range(i+1):
            print(i, j, i**2+j**2, np.sqrt(i**2+j**2))
