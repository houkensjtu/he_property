from math import sqrt, exp

class Helium(object):
    """
    The Helium-4 gas model. Most property equations come from NIST.
    """

    # The coef of 32 terms equation of state.
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

    def __init__(self, temperature=300.0, pressure=1.0e6):
        """ Initialize a Helium object. 
        Default temperature = 300 K
        Default pressure = 1.0e6 pa
        Default density = 1.0
        when no parameter is given.
        """
        self.temp = temperature
        self.p    = pressure
        self.pp   = self.p / 1.0e6
        
        self.rp   = 1.0
        self.h    = 0.0
        self.s    = 0.0


    def rpcalcul(self):
        """
        Return density(kg/m^3) based on temperature(K) ,pressure(MPa) and 
        an initial guess rp.
        """
        temp = self.temp
        pp   = self.pp
        rp   = self.rp
        
        _TAO = 0.0033033259
        _R   = 0.00831434
        _G   = self._G
        
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
        # Notice this is not a polynomial of rp, because there is a exp().
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

        self.rp = rp0 * 4.0026
        
        return rp0 * 4.0026


    def prop2(self,t,d):
        """
        Return SD based on temperature and density.
        """
        
        gamma = -.33033259e-02
        r     = 0.00831434
        f     = exp(gamma * d**2)

        G1 = f/(2.*gamma)
        G2 = (f*d**2 - 2.*G1)/(2. * gamma)
        G3 = (f*d**4 - 4.*G2)/(2. * gamma)
        G4 = (f*d**6 - 6.*G3)/(2. * gamma)
        G5 = (f*d**8 - 8.*G4)/(2. * gamma)
        G6 = (f*d**10 - 10.*G5)/(2. * gamma)

        X = [None] * 32
        X[0] = -d
        X[1] = -d / (2. * sqrt(t))
        X[2] = 0.
        X[3] = d / t**2
        X[4] = 2.*d / t**3
        X[5] = -d**2 / 2.
        X[6] = 0.
        X[7] = d**2 / (2.* t**2)
        X[8] = d**2 / t**3
        X[9] = -d**3 / 3.
        X[10] = 0.
        X[11] = d**3 / (3.* t**2)
        X[12] = 0.
        X[13] = d**5 / (5.* t**2)
        X[14] = 2.*d**5/(5.* t**3)
        X[15] = d**6 / (6.* t**2)
        X[16] = d**7 / (7.* t**2)
        X[17] = 2.*d**7 / (7. * t**3)
        X[18] = d**8 / (4.* t**3)
        X[19] = 2.*G1 / t**3
        X[20] = 3.*G1 / t**4
        X[21] = 2.*G2 / t**3
        X[22] = 4.*G2 / t**5
        X[23] = 2.*G3 / t**3
        X[24] = 3.*G3 / t**4
        X[25] = 2.*G4 / t**3
        X[26] = 4.*G4 / t**5
        X[27] = 2.*G5 / t**3
        X[28] = 3.*G5 / t**4
        X[29] = 2.*G6 / t**3
        X[30] = 3.*G6 / t**4
        X[31] = 4.*G6 / t**5
        
        sd = 0
        # from 0 -> 31
        for i in range(0,32):
            sd += self._G[i] * X[i]

        return sd


    def prop3(self,t,d):
        """
        Return UD based on temperature and density.
        """
        gamma = -.33033259e-02
        r     = 0.00831434
        f = exp(gamma * d**2)

        G1 = f / (2. * gamma)
        G2 = (f * d**2 - 2.*G1) / (2. * gamma)
        G3 = (f * d**4 - 4.*G2) / (2. * gamma)
        G4 = (f * d**6 - 6.*G3) / (2. * gamma)
        G5 = (f * d**8 - 8.*G4) / (2. * gamma)
        G6 = (f * d**10 - 10.*G5) / (2. * gamma)

        X = [None] * 32
        X[0] = d * t
        X[1] = d * sqrt(t)
        X[2] = d
        X[3] = d / t
        X[4] = d / t**2
        X[5] = d**2 * t / 2.
        X[6] = d**2 / 2.
        X[7] = d**2 / (2.*t)
        X[8] = d**2 / (2.*t**2)
        X[9] = d**3 * t / 3.
        X[10] = d**3 / 3.
        X[11] = d**3 / (3.*t)
        X[12] = d**4 / 4.
        X[13] = d**5 / (5.*t)
        X[14] = d**5 / (5.*t**2)
        X[15] = d**6 / (6.*t)
        X[16] = d**7 / (7.*t)
        X[17] = d**7 / (7.*t**2)
        X[18] = d**8 / (8.*t**2)
        X[19] = G1 / t**2
        X[20] = G1 / t**3
        X[21] = G2 / t**2
        X[22] = G2 / t**4
        X[23] = G3 / t**2
        X[24] = G3 / t**3
        X[25] = G4 / t**2
        X[26] = G4 / t**4
        X[27] = G5 / t**2
        X[28] = G5 / t**3
        X[29] = G6 / t**2
        X[30] = G6 / t**3
        X[31] = G6 / t**4
        
        ud = 0
        for i in range(0,32):
            ud += self._G[i] * X[i]
            
        return ud

    def cpi(self,t):
        """ Return specific heat capacity of helium (J/mol*K). 4-300K."""
        rr = 8.31434
        g  = [0.0, 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 2.5, 0.0, 0.0,]
        cpi = g[3] * t
        cpi *= rr
        return cpi

    def enthm(self):
        """ Return enthalpy(J/kg) based on pressure(Pa), density(kg/m^3) and temperature."""
        s0 = 1964.2
        h0 = 15273.0
        r  = 8.31434e-3

        dd1 = self.rp/4.00246
        tt  = self.temp
        pp  = self.pp
        
        sd  = self.prop2(tt,dd1)
        ud  = self.prop3(tt,dd1)
        
        dd0 = 0.0
        
        s0  = self.prop2(tt,dd0)
        u0  = self.prop3(tt,dd0)

        enthm = tt * (sd-s0)*1000 + (ud-u0)*1000 + self.cpi(tt) + (pp/dd1 - r * tt)*1.0e3
        enthm = enthm / 4.0026*1.0e3
        enthm = enthm + h0
        self.h = enthm
        return enthm

    
def main():
    p = 0.2e6
    for t in range(3, 300, 1):
        he = Helium(temperature=t, pressure=p)
        d  = he.rpcalcul()
        h  = he.enthm()
        print ("Temp:",t," Pressure:", p,"Density:", d,"Enthalpy:", h)

if __name__ == "__main__":
    main()
