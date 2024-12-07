#!/usr/bin/env python3
# Name: Justin Jang (jjang12)
# Group Members: Aster Lathbury(mlathbur), Kimberly Magpantay(klmagpan), Faiz Khan(faahkhan)
'''
Program docstring goes here

Read a set of 3 atomic coordinates, 2 and return the bond length between N-C, bond length between N-Ca, and the bond angle of C-N-Ca (think geometry bond angles from Chem1b)

For example, if I enter the following coordinates (notice.. they are all on one line !!!) :<br>

 <br>
C = (39.447, 94.657, 11.824) N = (39.292, 95.716, 11.027) Ca = (39.462, 97.101, 11.465)
then the program will output the following three lines:<br>
N-C bond length = 1.33<br>
N-Ca bond length = 1.46<br>
C-N-Ca bond angle = 124.0<br> 
'''

import math
import re
class Triad :
    """
    Calculate angles and distances among a triad of points.
 
    Author: David Bernick
    Date: March 21, 2013
    Points can be supplied in any dimensional space as long as they are consistent.
    Points are supplied as tupels in n-dimensions, and there should be three
    of those to make the triad. Each point is positionally named as p,q,r
    and the corresponding angles are then angleP, angleQ and angleR.
    Distances are given by dPQ(), dPR() and dQR()
 
    Required Modules: math , re
    initialized: 3 positional tuples representing Points in n-space
             p1 = Triad( p=(1,0,0), q=(0,0,0), r=(0,1,0) )
    attributes: p,q,r the 3 tuples representing points in N-space
    methods:  angleP(), angleR(), angleQ() angles measured in radians
          dPQ(), dPR(), dQR() distances in the same units of p,q,r
 
    """
 
    def __init__(self,p,q,r) :
        """ Construct a Triad. 
        
        Example construction:
            p1 = Triad( p=(1.,0.,0.), q=(0.,0.,0.), r=(0.,0.,0.) ). 
        """
        self.p = p
        self.q = q
        self.r = r
        
# private helper methods
    def d2 (self,a,b) : # calculate squared distance of point a to b
        return float(sum((ia-ib)*(ia-ib)  for  ia,ib in zip (a,b)))
    
    def dot (self,a,b) : # dotProd of standard vectors a,b
        return float(sum(ia*ib for ia,ib in zip(a,b)))
    
    def ndot (self,a,b,c) : # dotProd of vec. a,c standardized to b
        return float(sum((ia-ib)*(ic-ib) for ia,ib,ic in zip (a,b,c)))
    
# calculate lengths(distances) of segments PQ, PR and QR
    def dPQ (self):
        """ Provides the distance between point p and point q """
        return math.sqrt(self.d2(self.p,self.q))
    
    def dPR (self):
        """ Provides the distance between point p and point r """
        return math.sqrt(self.d2(self.p,self.r))
    
    def dQR (self):
        """ Provides the distance between point q and point r """
        return math.sqrt(self.d2(self.q,self.r))
    
    def angleP (self) :
        """ Provides the angle made at point p by segments pq and pr (radians). """
        return math.acos(self.ndot(self.q,self.p,self.r) / 
                         math.sqrt(self.d2(self.q,self.p) * self.d2(self.r,self.p)))
    
    def angleQ (self) :
        """ Provides the angle made at point q by segments qp and qr (radians). """
        return math.acos(self.ndot(self.p,self.q,self.r) /
                         math.sqrt(self.d2(self.p,self.q) * self.d2(self.r,self.q)))
 
    def angleR (self) :
        """ Provides the angle made at point r by segments rp and rq (radians). """
        return math.acos(self.ndot(self.p,self.r,self.q) /
                         math.sqrt(self.d2(self.p,self.r) * self.d2(self.q,self.r)))

#user-added code starts below:

def main():
    # Extract coordinates from input string
    coord = re.findall(r'[-+]?\d*\.\d+,\s*[-+]?\d*\.\d+,\s*[-+]?\d*\.\d+', input(''))
    #use of regex, same pattern is repeated 3 times because the input string has 3 positional tuples, will match a string that contains 3 comma-separated set of numbers.
    #findall() searches for all instances in the string

    #convert comma-separated set of numbers into tuples of floats
    C = tuple(map(float, coord[0].split(', ')))
    N = tuple(map(float, coord[1].split(', ')))
    Ca = tuple(map(float, coord[2].split(', ')))

    # Create Triad object
    triad = Triad(p=C, q=N, r=Ca)
    
    # Calculate bond lengths and angles from template code above
    NC_length = triad.dPQ()
    NCa_length = triad.dPR()
    CNCa_angleDEGREES = math.degrees(triad.angleQ())  #in degrees, per prompt

    # Print the results in 3 lines
    print(f"N-C bond length = {NC_length:.2f}")
    print(f"N-Ca bond length = {NCa_length:.2f}")
    print(f"C-N-Ca bond angle = {CNCa_angleDEGREES:.1f}")

main()
 