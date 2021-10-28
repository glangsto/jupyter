#python function to calculate the angular distance between to coordinate pairs
#
import math
import numpy as np

def distance(coord1, coord2):
    """
    angulardistance takes a pair of lon,lat coordinates and
    computes the angular distance between the coordinates.
    Inputs and outputs are in degrees
    """
    lat1 = float(coord1[1])
    long1 = float(coord1[0])
    lat2 = float(coord2[1])
    long2 = float(coord2[0])
    # Convert latitude and longitude to spherical coordinates in radians.
    degrees_to_radians = math.pi/180.0
    phi1 = (90.0 - lat1)*degrees_to_radians
    phi2 = (90.0 - lat2)*degrees_to_radians
    theta1 = long1*degrees_to_radians
    theta2 = long2*degrees_to_radians
    # Compute spherical distance from spherical coordinates.
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + 
           math.cos(phi1)*math.cos(phi2))
    # compute the arc-cosign to get angular distance in radians
    arc = math.acos( cos )
    # Convert back to degrees
    distanceDeg = arc / degrees_to_radians
    return distanceDeg
