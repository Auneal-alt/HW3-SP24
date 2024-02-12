# Note code from Smays help video to plug into HW3b.py

def Simpson(fcn, args, a, b, npoints=21):
    """Simposn integraition code for integrating or finding the area below a curve"""

    area = 0  # intital value for the integral
    m = npoints
    n = 2 * m  # THIS ensures and even number of panels
    xL = min(a, b)
    xR = max(a, b)

    if xL == xR:
        return 0

    h = (xR - xL) / n
    x = xL
    area = (fcn(xL, args) + fcn(xR, args))

    for j in range(1, n):  # counts from 1 to 2*m-1
        x = j * h + xL  # updates the x position for evaluating the function
        if not j % 2 == 0:  # the odds
            area += 4 * fcn((x, args))
        else:  # the evens
            area += 2 * fcn((x, args))
    return (h / 3.0) * area  # finally, return the value for the integral
