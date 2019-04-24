import numpy as np
import matplotlib.pyplot as plt

def main():
    NACA4digit("2412")


def NACA4digit(airfoil_number, points=1000, closed_trailing_edge=True, save_to_file=False):
    """Compute and plot coordinates of a NACA 4-digit airfoil.
    
    Arguments:
        number: Name of the requested 4-digt NACA airfoil entered as a string,
                i.e. '2412'
    
    Optional arguments:
        points: Number of desired airfoil coordinates
                (half(N) for upper and half(N) for lower surface)
        closed_trailing_edge: The trailing edge has a small but finite
                            thickness when using equations in [1] unaltered.
        save_to_file: Saves the coordinates in a CSV file.
    
    Returns:
        A matrix with two columns pertaining to the x and y-coordinates,
        respectively. The sequence of coordinates is clockwise, starting at the
        trailing edge.
    
    Raises:
        ValueError: Airfoil number does not have four digits or N is negative
        
    References:
        https://en.wikipedia.org/wiki/NACA_airfoil
    
    """

    ## Check whether input parameters are valid or not
    if not (0 < int(airfoil_number) < 10000 and points > 0):
        raise ValueError("Invalid input.")

    c = 1       # by default chord(c) is set to 1; so x = c, x/c = 1

    ## 4-digit airfoil properties
    M = int(airfoil_number[0]) * (c/100)    # max camber
    P = int(airfoil_number[1]) * (c/10)     # max camber location from leading edge
    T = int(airfoil_number[2:4]) * (c/100)  # max thickness

    ## coordinates points
    N = points//2

    ## spacing of coordinates along the x-axis
    x = np.linspace(0, c, N)

    ## The trailing edge has a small but finite thickness by default. The gap
    # can be closed by utilizing a slightly different equation.
    if closed_trailing_edge:
        a4 = -0.1036   # Closed trailing edge
    else:
        a4 = -0.1015   # Open trailing edge

    ## constants
    a0 = 0.2969
    a1 = -0.1260
    a2 = -0.3516
    a3 = 0.2843

    ## Computing the y-coordinates of the camber line.
    # Camber & Gradient
    fwd_x = x[x < P]
    aft_x = x[x >= P]
    # camber calcultation
    if 0 < P < c:
        fwd_camber_yc = M / P**2 * (2 * P * fwd_x - np.power(fwd_x, 2))
        fwd_camber_dyc_dx = 2*M / P**2 * (P - fwd_x)
        aft_camber_yc = M / (1 - P)**2 * ((1 - 2 * P) + 2 * P * aft_x - np.power(aft_x, 2))
        aft_camber_dyc_dx = 2*M / (1-P)**2 * (P - aft_x)
        camber_yc = np.append(fwd_camber_yc, aft_camber_yc)   # complete mean camber line
        dyc_dx_camber = np.append(fwd_camber_dyc_dx, aft_camber_dyc_dx)
    else:
        camber_yc = np.zeros(np.size(x))   # cases where max_camber is located on leading or triling edge
        dyc_dx_camber = np.zeros(np.size(x))
    #gradient calculation
    theta = np.arctan(dyc_dx_camber)
    
    ## Thickness Distribution
    yt = 5*T*((a0*np.sqrt(x)) + (a1*x) + (a2*(x**2)) + (a3*(x**3)) + (a4*(x**4)))

    ## Upper surface points
    xu = x - (yt * np.sin(theta))
    yu = camber_yc + (yt* np.cos(theta))

    ## Lower surface points
    xl = x + (yt * np.sin(theta))
    yl = camber_yc - (yt * np.cos(theta))
    
    ## coordinates
    x_cord = np.append(xu, xl)
    y_cord = np.append(yu, yl)
    coordinates = np.column_stack((x_cord, y_cord))
    
    ## Saving data to CSV file
    if save_to_file:
        np.savetxt(f"NACA{airfoil_number}_plotting_data.csv", coordinates, delimiter=',')

    ## Plot the airfoil
    plt.plot(x_cord, y_cord)
    plt.grid()
    plt.axis('equal')
    return plt.show()


if __name__ == '__main__':
    main()

