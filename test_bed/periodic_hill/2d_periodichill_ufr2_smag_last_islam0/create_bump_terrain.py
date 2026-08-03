import numpy as np
from scipy.interpolate import PPoly
import matplotlib.pyplot as plt
from bisect import bisect_right
import pandas as pd


def calcz(x,xb,coeffs):
    ind = bisect_right(xb,x)-1
    c = coeffs[:,ind]
    return c[3] + c[2]*x + c[1]*x*x + c[0]*x*x*x

def hill(x):
    h = 28.0
    xstar = x * h
    xstar[xstar > 128] = 252 - xstar[xstar > 128]
    ystar = np.zeros(x.shape)
    idx = (0.0 <= xstar) & (xstar < 9.0)
    ystar[idx] = np.minimum(
        28 * np.ones(x[idx].shape),
        2.800000000000e01
        + 0.000000000000e00 * xstar[idx]
        + 6.775070969851e-03 * xstar[idx] ** 2
        - 2.124527775800e-03 * xstar[idx] ** 3,
    )
    idx = (9.0 <= xstar) & (xstar < 14.0)
    ystar[idx] = (
        2.507355893131e01
        + 9.754803562315e-01 * xstar[idx]
        - 1.016116352781e-01 * xstar[idx] ** 2
        + 1.889794677828e-03 * xstar[idx] ** 3
    )
    idx = (14.0 <= xstar) & (xstar < 20.0)
    ystar[idx] = (
        2.579601052357e01
        + 8.206693007457e-01 * xstar[idx]
        - 9.055370274339e-02 * xstar[idx] ** 2
        + 1.626510569859e-03 * xstar[idx] ** 3
    )
    idx = (20.0 <= xstar) & (xstar < 30.0)
    ystar[idx] = (
        4.046435022819e01
        - 1.379581654948e00 * xstar[idx]
        + 1.945884504128e-02 * xstar[idx] ** 2
        - 2.070318932190e-04 * xstar[idx] ** 3
    )
    idx = (30.0 <= xstar) & (xstar < 40.0)
    ystar[idx] = (
        1.792461334664e01
        + 8.743920332081e-01 * xstar[idx]
        - 5.567361123058e-02 * xstar[idx] ** 2
        + 6.277731764683e-04 * xstar[idx] ** 3
    )
    idx = (40.0 <= xstar) & (xstar < 50.0)
    ystar[idx] = np.maximum(
        np.zeros(x[idx].shape),
        5.639011190988e01
        - 2.010520359035e00 * xstar[idx]
        + 1.644919857549e-02 * xstar[idx] ** 2
        + 2.674976141766e-05 * xstar[idx] ** 3,
    )

    return ystar / h


def main():

    xplusmax = 9.0
    yplusmax = 4.5
    h = 28.0
    hm = h/1000.0
    xmax = xplusmax*h
    ymax = yplusmax*h
    zadd = 0.00096591 # One cell height where H_z=0.085 and N_z=88

    ngridx = 100
    ngridy = 50

    xnew = np.linspace(0, xplusmax, ngridx)
    znew = np.array(hill(xnew))
    ynew = np.linspace(0,yplusmax,ngridy)


    # Plot to verify
    plt.figure(figsize=(8, 5))
    plt.plot(xnew, znew, color='blue', label='Cubic Poly Spline')

    # Add formatting
    plt.title('Hill')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.legend()
    plt.axis('equal') 
    plt.grid(True)
    plt.savefig('hill.png')

    # Output csv

    outdf = pd.DataFrame({
        'x':xnew*hm,
        'z':znew*hm
    })

    outdf.to_csv('periodic_hill.csv', index=False)

    # Output kynema file
    zall = np.array([])
    for x in range(xnew.size):
        zall = np.append(zall,np.ones(ngridy)*znew[x])

    output_array = np.array([])
    output_array = np.append(output_array,xnew*hm)
    output_array = np.append(output_array,ynew*hm)
    output_array = np.append(output_array,zall*hm + zadd)

    np.savetxt("terrain.txt", [xnew.size,ynew.size], fmt="%d")
    with open("terrain.txt", "a") as f:
        np.savetxt(f, output_array , fmt="%f")


if __name__ == "__main__":
    main()

