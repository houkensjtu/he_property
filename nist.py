class Helium(object):
    """
    The Helium-4 gas model. Most property equations come from NIST.
    """    
    def __init__(self, temperature = 300.0, pressure = 1.0e6, density = 1.0):
        """ Initialize a Helium object. 
        Default temperature = 300 K
        Default pressure = 1.0e6 pa
        Default density = 1.0
        when no parameter is given.
        """
        self.temp = temperature
        self.p    = pressure
        self.pp   = self.p / 1.0e6
        self.rp   = density

    def rpcalcul(self):
        """
        Return density(kg/m^3) based on temperature(K) ,pressure(Pa) and 
        an initial guess rp.
        """
        from math import sqrt, exp
        temp = self.temp
        pp   = self.pp
        rp   = self.rp
        
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

        # f is the function of pressure, density and temperature.
        # when f(rp,pp,t)=0
        # means rp is the density solution for pp and t    
        f = lambda rp : pp - ( rp*a1+rp*rp*a2+rp**3*a3+rp**4*a4+rp**5*a5+rp**6*a6
                             +rp**7*a7+rp**8*a8+rp**9*a9
                             +(rp**3*b3+rp**5*b5+rp**7*b7+rp**9*b9+rp**11*b11
                             +rp**13*b13) * exp(-_TAO*rp*rp))

        # fp is the derivative df/d(rp) of f
        fp = lambda rp : -(a1+2.0e0*rp*a2+3.0e0*rp**2*a3+4.0e0*rp**3*a4
                         +5.0e0*rp**4*a5+6.0e0*rp**5*a6+7.0e0*rp**6*a7
                         +8.0e0*rp**7*a8+9.0e0*rp**8*a9+(3.0e0*rp**2*b3
                         +5.0e0*rp**4*b5+7.0e0*rp**6*b7+9.0e0*rp**8*b9
                         +11.0e0*rp**10*b11+13.0e0*rp**12*b13) * exp(-_TAO*rp*rp)
                         +(rp**3*b3+rp**5*b5+rp**7*b7+rp**9*b9+rp**11*b11
                         +rp**13*b13) * exp(-_TAO*rp*rp)*(-2.0e0*_TAO*rp))

        # Give an initial guess of rp based on temperature.
        rp0 = 0.0
        if self.temp > 8.0:
            rp0 = rp / 4.0026
        else:
            rp0 = 160.0 / 4.0026

        # The Newton root-finding algorithm
        while abs(f(rp0)) > EPS:
            f0 = f(rp0)
            fp0 = fp(rp0)
            rp0 = rp0 - f0 / fp0
            if rp0 <= 0.0:
                rp0 = rp0 + 0.1 * rp0

        return rp0 * 4.0026


def main():
    for t in range(4, 300, 10):
        he = Helium(temperature=t, pressure=.5e6)
        print (t, he.rpcalcul())

if __name__ == "__main__":
    main()
