"""Response and frequency data with optional Matplotlib rendering."""
import numpy as np
import symox as S
from .lti import _tf


def _numeric(system):
    system=_tf(system)
    if system.to_expr().free_symbols-{system.var}:raise ValueError('numerical data require numeric system parameters')
    return system


def _poly_roots(poly):
    def _eval(values):return [float(v) if v.is_real else complex(v) for v in values]
    if poly.domain in (S.QQ,S.ZZ):return _eval(poly.all_roots())
    return _eval(poly.nroots())


def pole_zero_numerical_data(system):
    """Return (zeros, poles) as Python floats/complexes, like the reference."""
    system=_numeric(system)
    return _poly_roots(S.Poly(system.num,system.var)),_poly_roots(S.Poly(system.den,system.var))


def _response(system,order,prec=8,lower_limit=0,upper_limit=10,slope=1,**kwargs):
    system=_numeric(system);count=int(kwargs.pop('n',1000));kwargs.pop('adaptive',None)
    if kwargs:raise TypeError('unknown response options: '+', '.join(kwargs))
    if count<2 or upper_limit<=lower_limit or lower_limit<0:raise ValueError('time bounds must be ordered and nonnegative, with at least two samples')
    if order==0 and not system.is_strictly_proper:raise ValueError('an impulse feedthrough contains a Dirac delta; finite numerical samples require a strictly proper transfer function')
    t=S.Dummy('response_time',real=True)
    expression=S.inverse_laplace_transform(system.to_expr()/system.var**order,system.var,t).subs(S.Heaviside(t),1)
    if order==2:expression*=S.sympify(slope)
    times=np.linspace(float(lower_limit),float(upper_limit),count)
    try:values=np.array([float(S.re(expression.subs(t,S.Float(float(value))).evalf(max(15,prec)))) for value in times])
    except (TypeError,ValueError) as error:raise NotImplementedError('response could not be evaluated from its inverse Laplace transform') from error
    return times,values


def step_response_numerical_data(system,prec=8,lower_limit=0,upper_limit=10,**kwargs):return _response(system,1,prec,lower_limit,upper_limit,**kwargs)
def impulse_response_numerical_data(system,prec=8,lower_limit=0,upper_limit=10,**kwargs):return _response(system,0,prec,lower_limit,upper_limit,**kwargs)
def ramp_response_numerical_data(system,slope=1,prec=8,lower_limit=0,upper_limit=10,**kwargs):return _response(system,2,prec,lower_limit,upper_limit,slope,**kwargs)


def _frequency_data(system,initial_exp,final_exp,freq_unit,**kwargs):
    system=_numeric(system);count=int(kwargs.pop('n',1000));kwargs.pop('adaptive',None)
    if kwargs:raise TypeError('unknown frequency options: '+', '.join(kwargs))
    if freq_unit not in ('rad/sec','Hz'):raise ValueError('frequency unit must be rad/sec or Hz')
    frequencies=np.logspace(float(initial_exp),float(final_exp),count)
    omega=frequencies*(2*np.pi if freq_unit=='Hz' else 1)
    values=np.array([complex(system.eval_frequency(S.I*S.Float(float(w))).evalf()) for w in omega])
    return frequencies,values


def bode_magnitude_numerical_data(system,initial_exp=-5,final_exp=5,freq_unit='rad/sec',**kwargs):
    f,v=_frequency_data(system,initial_exp,final_exp,freq_unit,**kwargs)
    return f,20*np.log10(np.abs(v))


def bode_phase_numerical_data(system,initial_exp=-5,final_exp=5,freq_unit='rad/sec',phase_unit='rad',phase_unwrap=True,**kwargs):
    if phase_unit not in ('rad','deg'):raise ValueError('phase unit must be rad or deg')
    f,v=_frequency_data(system,initial_exp,final_exp,freq_unit,**kwargs);phase=np.angle(v)
    if phase_unwrap:phase=np.unwrap(phase)
    return f,phase*(180/np.pi if phase_unit=='deg' else 1)


def _plot(data,xlabel,ylabel,show=True,grid=True,color='b',logx=False,**kwargs):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots();(ax.semilogx if logx else ax.plot)(*data,color=color,**kwargs);ax.set_xlabel(xlabel);ax.set_ylabel(ylabel);ax.grid(grid)
    if show:plt.show()
    return fig


def pole_zero_plot(system,pole_color='blue',pole_markersize=10,zero_color='orange',zero_markersize=7,grid=True,show_axes=True,show=True,**kwargs):
    import matplotlib.pyplot as plt
    zeros,poles=pole_zero_numerical_data(system);fig,ax=plt.subplots()
    ax.scatter([v.real for v in poles],[v.imag for v in poles],color=pole_color,marker='x',s=pole_markersize**2,label='poles')
    ax.scatter([v.real for v in zeros],[v.imag for v in zeros],edgecolors=zero_color,facecolors='none',s=zero_markersize**2,label='zeros')
    ax.grid(grid);ax.set_xlabel('Real');ax.set_ylabel('Imaginary')
    if show_axes:ax.axhline(0,color='black',linewidth=0.5);ax.axvline(0,color='black',linewidth=0.5)
    if show:plt.show()
    return fig


def step_response_plot(system,color='b',prec=8,lower_limit=0,upper_limit=10,show_axes=False,grid=True,show=True,**kwargs):return _plot(step_response_numerical_data(system,prec,lower_limit,upper_limit,**kwargs),'Time','Step response',show,grid,color)
def impulse_response_plot(system,color='b',prec=8,lower_limit=0,upper_limit=10,show_axes=False,grid=True,show=True,**kwargs):return _plot(impulse_response_numerical_data(system,prec,lower_limit,upper_limit,**kwargs),'Time','Impulse response',show,grid,color)
def ramp_response_plot(system,slope=1,color='b',prec=8,lower_limit=0,upper_limit=10,show_axes=False,grid=True,show=True,**kwargs):return _plot(ramp_response_numerical_data(system,slope,prec,lower_limit,upper_limit,**kwargs),'Time','Ramp response',show,grid,color)
def bode_magnitude_plot(system,initial_exp=-5,final_exp=5,color='b',show_axes=False,grid=True,show=True,freq_unit='rad/sec',**kwargs):return _plot(bode_magnitude_numerical_data(system,initial_exp,final_exp,freq_unit,**kwargs),'Frequency ('+freq_unit+')','Magnitude (dB)',show,grid,color,True)
def bode_phase_plot(system,initial_exp=-5,final_exp=5,color='b',show_axes=False,grid=True,show=True,freq_unit='rad/sec',phase_unit='rad',phase_unwrap=True,**kwargs):return _plot(bode_phase_numerical_data(system,initial_exp,final_exp,freq_unit,phase_unit,phase_unwrap,**kwargs),'Frequency ('+freq_unit+')','Phase ('+phase_unit+')',show,grid,color,True)


def bode_plot(system,initial_exp=-5,final_exp=5,grid=True,show_axes=False,show=True,freq_unit='rad/sec',phase_unit='rad',phase_unwrap=True,**kwargs):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,1,sharex=True)
    axes[0].semilogx(*bode_magnitude_numerical_data(system,initial_exp,final_exp,freq_unit,**kwargs))
    axes[1].semilogx(*bode_phase_numerical_data(system,initial_exp,final_exp,freq_unit,phase_unit,phase_unwrap,**kwargs))
    axes[0].set_ylabel('Magnitude (dB)');axes[1].set_ylabel('Phase ('+phase_unit+')');axes[1].set_xlabel('Frequency ('+freq_unit+')')
    for ax in axes:ax.grid(grid)
    if show:plt.show()
    return fig


def nyquist_plot_expr(system):
    system=_tf(system);omega=S.Symbol('omega',real=True)
    real,imag=system.eval_frequency(S.I*omega).as_real_imag()
    return S.simplify(real),S.simplify(imag),omega


def nichols_plot_expr(system):
    system=_tf(system);omega=S.Symbol('omega',positive=True);expression=system.eval_frequency(S.I*omega)
    return S.arg(expression)*180/S.pi,20*S.log(S.Abs(expression),10),omega


def nyquist_plot(system,initial_omega=0.01,final_omega=100,show=True,color='b',**kwargs):
    frequencies,values=_frequency_data(system,np.log10(initial_omega),np.log10(final_omega),'rad/sec',**kwargs)
    return _plot((values.real,values.imag),'Real','Imaginary',show,True,color)


def nichols_plot(system,initial_omega=0.01,final_omega=100,show=True,color='b',**kwargs):
    frequencies,values=_frequency_data(system,np.log10(initial_omega),np.log10(final_omega),'rad/sec',**kwargs)
    return _plot((np.unwrap(np.angle(values))*180/np.pi,20*np.log10(np.abs(values))),'Phase (deg)','Magnitude (dB)',show,True,color)


__all__=['pole_zero_numerical_data','pole_zero_plot','step_response_numerical_data','step_response_plot','impulse_response_numerical_data','impulse_response_plot','ramp_response_numerical_data','ramp_response_plot','bode_magnitude_numerical_data','bode_phase_numerical_data','bode_magnitude_plot','bode_phase_plot','bode_plot','nyquist_plot_expr','nyquist_plot','nichols_plot_expr','nichols_plot']
