"""Monochromatic transverse waves and coherent superposition."""
import symox as S
from ..units import speed_of_light


class TWave:
    def __init__(self,amplitude,frequency=None,phase=0,time_period=None,n=None):
        if frequency is None and time_period is None:raise ValueError('frequency or time period is required')
        self.amplitude,self.phase=S.sympify(amplitude),S.sympify(phase)
        self.frequency=S.sympify(frequency) if frequency is not None else 1/S.sympify(time_period)
        if time_period is not None and S.simplify(self.frequency*S.sympify(time_period)-1).is_zero is False:raise ValueError('frequency and period disagree')
        self.n=S.Symbol('n') if n is None else S.sympify(n)
    @property
    def time_period(self):return 1/self.frequency
    @property
    def angular_velocity(self):return 2*S.pi*self.frequency
    @property
    def speed(self):return S.sympify(speed_of_light)/self.n
    @property
    def wavelength(self):return self.speed/self.frequency
    @property
    def wavenumber(self):return 2*S.pi/self.wavelength
    def equation(self,form='cos',x=None,t=None):
        x=S.Symbol('x') if x is None else S.sympify(x);t=S.Symbol('t') if t is None else S.sympify(t)
        phase=self.wavenumber*x-self.angular_velocity*t+self.phase
        if form=='cos':return self.amplitude*S.cos(phase)
        if form=='sin':return self.amplitude*S.sin(phase)
        if form=='exp':return self.amplitude*S.exp(S.I*phase)
        raise ValueError('wave form must be cos, sin or exp')
    def __add__(self,other):
        if not isinstance(other,TWave) or self.frequency!=other.frequency or self.n!=other.n:raise ValueError('superposition requires equal frequencies and refractive indices')
        real=self.amplitude*S.cos(self.phase)+other.amplitude*S.cos(other.phase)
        imag=self.amplitude*S.sin(self.phase)+other.amplitude*S.sin(other.phase)
        return TWave(S.sqrt(real*real+imag*imag),self.frequency,S.atan2(imag,real),n=self.n)


__all__=['TWave']
