"""Jones and Stokes polarization, retarders, filters and beam splitters."""
import symox as S


def jones_vector(psi,chi):
    psi,chi=map(S.sympify,(psi,chi));return S.Matrix([S.cos(psi)*S.cos(chi)-S.I*S.sin(psi)*S.sin(chi),S.sin(psi)*S.cos(chi)+S.I*S.cos(psi)*S.sin(chi)])


def stokes_vector(psi,chi,p=1,I=1):
    psi,chi,p,I=map(S.sympify,(psi,chi,p,I));return I*S.Matrix([1,p*S.cos(2*chi)*S.cos(2*psi),p*S.cos(2*chi)*S.sin(2*psi),p*S.sin(2*chi)])


def jones_2_stokes(e):
    e=S.Matrix(e)
    if e.shape!=(2,1):raise ValueError('Jones vector must have two components')
    # Reference formula, without extra simplification.
    ex,ey=e
    return S.Matrix([S.Abs(ex)**2+S.Abs(ey)**2,S.Abs(ex)**2-S.Abs(ey)**2,
                     2*S.re(ex*ey.conjugate()),-2*S.im(ex*ey.conjugate())])


def linear_polarizer(theta=0):
    theta=S.sympify(theta);c,s=S.cos(theta),S.sin(theta);return S.Matrix([[c*c,c*s],[c*s,s*s]])


def phase_retarder(theta=0,delta=0):
    # Reference matrix and phase factor, entry for entry.
    theta,delta=map(S.sympify,(theta,delta))
    R=S.Matrix([[S.cos(theta)**2+S.exp(S.I*delta)*S.sin(theta)**2,(1-S.exp(S.I*delta))*S.cos(theta)*S.sin(theta)],
                [(1-S.exp(S.I*delta))*S.cos(theta)*S.sin(theta),S.sin(theta)**2+S.exp(S.I*delta)*S.cos(theta)**2]])
    return R*S.exp(-S.I*delta/2)


def half_wave_retarder(theta):return phase_retarder(theta,S.pi)
def quarter_wave_retarder(theta):return phase_retarder(theta,S.pi/2)
def transmissive_filter(T):return S.sqrt(S.sympify(T))*S.eye(2)
def reflective_filter(R):return S.sqrt(S.sympify(R))*S.diag(1,-1)


def mueller_matrix(J):
    J=S.Matrix(J)
    if J.shape!=(2,2):raise ValueError('Jones matrix must be 2 by 2')
    # Reference: simplify(A*TensorProduct(J, J.conjugate())*A.inv()).
    A=S.Matrix([[1,0,0,1],[1,0,0,-1],[0,1,1,0],[0,-S.I,S.I,0]])
    return S.simplify(A*S.kronecker_product(J,J.conjugate())*A.inv())


def polarizing_beam_splitter(Tp=1,Rs=1,Ts=0,Rp=0,phia=0,phib=0):
    Tp,Rs,Ts,Rp,phia,phib=map(S.sympify,(Tp,Rs,Ts,Rp,phia,phib))
    return S.Matrix([[S.sqrt(Tp),0,S.I*S.sqrt(Rp),0],[0,S.sqrt(Ts),0,-S.I*S.sqrt(Rs)*S.exp(S.I*phia)],[S.I*S.sqrt(Rp),0,S.sqrt(Tp),0],[0,-S.I*S.sqrt(Rs)*S.exp(S.I*phib),0,S.sqrt(Ts)]])


__all__=['jones_vector','stokes_vector','jones_2_stokes','linear_polarizer','phase_retarder','half_wave_retarder','quarter_wave_retarder','transmissive_filter','reflective_filter','mueller_matrix','polarizing_beam_splitter']
