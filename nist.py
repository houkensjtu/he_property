class Helium(object):
    """
    The Helium-4 gas model. Most property equations come from NIST.
    """    
    def __init__(self, temperature, pressure):
        self.temp = temperature
        self.p    = pressure
        self.pp   = self.p / 1.0e6

    def density(self):
        
        from math import sqrt
        temp = self.temp
        pp   = self.pp
        
        _TAO = 0.0033033259
        _R   = 0.00831434
        _G   = [  .4558980227431e-04,  .1260692007853e-02,
             -.7139657549318e-02,  .9728903861441e-02, -.1589302471562e-01,
             .1454229259623e-05, -.4708238429298e-04,  .1132915223587e-02,
             .2410763742104e-02, -.5093547838381e-08,  .2699726927900e-05,
             -.3954146691114e-04,  .1551961438127e-08,  .1050712335785e-07,
             -.5501158366750e-07, -.1037673478521e-09,  .6446881346448e-12,
             .3298960057071e-10, -.3555585738784e-12, -.6885401367690e-02,
             .9166109232806e-02, -.6544314242937e-05, -.3315398880031e-04,
             -.2067693644676e-07,  .3850153114958e-07, -.1399040626999e-10,
             -.1888462892389e-11, -.4595138561035e-14,  .6872567403738e-14,
             -.6097223119177e-18, -.7636186157005e-17,  .3848665703556e-17,
        ]
        
        a1 = _R*temp
        a2 = _G[0]*temp+_G[1]*sqrt(temp)+_G[2]+_G[3]/temp+_G[4]/temp**2
        a3 = _G[5]*temp+_G[6]+_G[7]/temp+_G[8]/temp**2
        a4 = _G[9]*temp+_G[10]+_G[11]/temp
        a5 = _G[12]
        a6 = _G[13]/temp+_G[14]/temp**2
        a7 = _G[15]/temp
        a8 = _G[16]/temp+_G[17]/temp**2
        a9 = _G[18]/temp**2
        b3  = _G[19]/temp**2+_G[20]/temp**3
        b5  = _G[21]/temp**2+_G[22]/temp**4
        b7  = _G[23]/temp**2+_G[24]/temp**3
        b9  = _G[25]/temp**2+_G[26]/temp**4
        b11 = _G[27]/temp**2+_G[28]/temp**3
        b13 = _G[29]/temp**2+_G[29]/temp**3+_G[31]/temp**4
        EPS = 1.0e-6
        
        f = lambda rp : pp - ( rp*a1+rp*rp*a2+rp**3*a3+rp**4*a4+rp**5*a5+rp**6*a6
                             +rp**7*a7+rp**8*a8+rp**9*a9
                             +(rp**3*b3+rp**5*b5+rp**7*b7+rp**9*b9+rp**11*b11
                             +rp**13*b13) ** (-_TAO*rp*rp))
        fp = lambda rp : -(a1+2.0e0*rp*a2+3.0e0*rp**2*a3+4.0e0*rp**3*a4
                         +5.0e0*rp**4*a5+6.0e0*rp**5*a6+7.0e0*rp**6*a7
                         +8.0e0*rp**7*a8+9.0e0*rp**8*a9+(3.0e0*rp**2*b3
                         +5.0e0*rp**4*b5+7.0e0*rp**6*b7+9.0e0*rp**8*b9
                         +11.0e0*rp**10*b11+13.0e0*rp**12*b13) ** (-_TAO*rp*rp)
                         +(rp**3*b3+rp**5*b5+rp**7*b7+rp**9*b9+rp**11*b11
                         +rp**13*b13) ** (-_TAO*rp*rp)*(-2.0e0*_TAO*rp))

        d = self.temp + self.p
        return d

def main():
    he = Helium(293.0, 2.50)
    print(he.density())

if __name__ == "__main__":
    main()
