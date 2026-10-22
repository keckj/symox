"""Snell refraction, Fresnel amplitudes and paraxial imaging relations.

The functions follow the reference implementation step by step, so that the
accepted inputs, total internal reflection branches, error messages and
returned types are the same.
"""
from symox import Float, I, oo, pi, zoo, S, Symbol, sympify, sqrt, acos, asin, atan2, cos, sin, tan, Matrix, cancel, Limit
from symox.geometry import Ray3D, Plane, intersection
from symox.utilities.iterables import is_sequence
from .medium import Medium

__all__ = ['refraction_angle', 'deviation', 'fresnel_coefficients', 'brewster_angle', 'critical_angle', 'lens_makers_formula', 'mirror_formula', 'lens_formula', 'hyperfocal_distance', 'transverse_magnification']

def refractive_index_of_medium(medium):
    if isinstance(medium, Medium):
        n = medium.refractive_index
    else:
        n = sympify(medium)
    return n

def refraction_angle(incident, medium1, medium2, normal=None, plane=None):
    n1 = refractive_index_of_medium(medium1)
    n2 = refractive_index_of_medium(medium2)
    try:
        angle_of_incidence = float(incident)
    except TypeError:
        angle_of_incidence = None
    try:
        critical_angle_ = critical_angle(medium1, medium2)
    except (ValueError, TypeError):
        critical_angle_ = None
    if angle_of_incidence is not None:
        if normal is not None or plane is not None:
            raise ValueError('Normal/plane not allowed if incident is an angle')
        if not 0.0 <= angle_of_incidence < pi * 0.5:
            raise ValueError('Angle of incidence not in range [0:pi/2)')
        if critical_angle_ and angle_of_incidence > critical_angle_:
            raise ValueError('Ray undergoes total internal reflection')
        return asin(n1 * sin(angle_of_incidence) / n2)
    return_ray = False
    if plane is not None and normal is not None:
        raise ValueError('Either plane or normal is acceptable.')
    if not isinstance(incident, Matrix):
        if is_sequence(incident):
            _incident = Matrix(incident)
        elif isinstance(incident, Ray3D):
            _incident = Matrix(incident.direction_ratio)
        else:
            raise TypeError('incident should be a Matrix, Ray3D, or sequence')
    else:
        _incident = incident
    if plane is not None:
        if not isinstance(plane, Plane):
            raise TypeError('plane should be an instance of geometry.plane.Plane')
        if isinstance(incident, Ray3D):
            return_ray = True
            intersection_pt = plane.intersection(incident)[0]
        _normal = Matrix(plane.normal_vector)
    elif not isinstance(normal, Matrix):
        if is_sequence(normal):
            _normal = Matrix(normal)
        elif isinstance(normal, Ray3D):
            _normal = Matrix(normal.direction_ratio)
            if isinstance(incident, Ray3D):
                intersection_pt = intersection(incident, normal)
                if len(intersection_pt) == 0:
                    raise ValueError("Normal isn't concurrent with the incident ray.")
                else:
                    return_ray = True
                    intersection_pt = intersection_pt[0]
        else:
            raise TypeError('Normal should be a Matrix, Ray3D, or sequence')
    else:
        _normal = normal
    eta = n1 / n2
    mag_incident = sqrt(sum((i ** 2 for i in _incident)))
    mag_normal = sqrt(sum((i ** 2 for i in _normal)))
    _incident /= mag_incident
    _normal /= mag_normal
    c1 = -_incident.dot(_normal)
    cs2 = 1 - eta ** 2 * (1 - c1 ** 2)
    if cs2.is_negative:
        return S.Zero
    drs = eta * _incident + (eta * c1 - sqrt(cs2)) * _normal
    drs = drs * mag_incident
    if not return_ray:
        return drs
    else:
        return Ray3D(intersection_pt, direction_ratio=drs)

def fresnel_coefficients(angle_of_incidence, medium1, medium2):
    if not 0 <= 2 * angle_of_incidence < pi:
        raise ValueError('Angle of incidence not in range [0:pi/2)')
    n1 = refractive_index_of_medium(medium1)
    n2 = refractive_index_of_medium(medium2)
    angle_of_refraction = asin(n1 * sin(angle_of_incidence) / n2)
    try:
        angle_of_total_internal_reflection_onset = critical_angle(n1, n2)
    except ValueError:
        angle_of_total_internal_reflection_onset = None
    if angle_of_total_internal_reflection_onset is None or angle_of_total_internal_reflection_onset > angle_of_incidence:
        R_s = -sin(angle_of_incidence - angle_of_refraction) / sin(angle_of_incidence + angle_of_refraction)
        R_p = tan(angle_of_incidence - angle_of_refraction) / tan(angle_of_incidence + angle_of_refraction)
        T_s = 2 * sin(angle_of_refraction) * cos(angle_of_incidence) / sin(angle_of_incidence + angle_of_refraction)
        T_p = 2 * sin(angle_of_refraction) * cos(angle_of_incidence) / (sin(angle_of_incidence + angle_of_refraction) * cos(angle_of_incidence - angle_of_refraction))
        return [R_p, R_s, T_p, T_s]
    else:
        n = n2 / n1
        R_s = cancel((cos(angle_of_incidence) - I * sqrt(sin(angle_of_incidence) ** 2 - n ** 2)) / (cos(angle_of_incidence) + I * sqrt(sin(angle_of_incidence) ** 2 - n ** 2)))
        R_p = cancel((n ** 2 * cos(angle_of_incidence) - I * sqrt(sin(angle_of_incidence) ** 2 - n ** 2)) / (n ** 2 * cos(angle_of_incidence) + I * sqrt(sin(angle_of_incidence) ** 2 - n ** 2)))
        return [R_p, R_s]

def deviation(incident, medium1, medium2, normal=None, plane=None):
    refracted = refraction_angle(incident, medium1, medium2, normal=normal, plane=plane)
    try:
        angle_of_incidence = Float(incident)
    except TypeError:
        angle_of_incidence = None
    if angle_of_incidence is not None:
        return float(refracted) - angle_of_incidence
    if refracted != 0:
        if isinstance(refracted, Ray3D):
            refracted = Matrix(refracted.direction_ratio)
        if not isinstance(incident, Matrix):
            if is_sequence(incident):
                _incident = Matrix(incident)
            elif isinstance(incident, Ray3D):
                _incident = Matrix(incident.direction_ratio)
            else:
                raise TypeError('incident should be a Matrix, Ray3D, or sequence')
        else:
            _incident = incident
        if plane is None:
            if not isinstance(normal, Matrix):
                if is_sequence(normal):
                    _normal = Matrix(normal)
                elif isinstance(normal, Ray3D):
                    _normal = Matrix(normal.direction_ratio)
                else:
                    raise TypeError('normal should be a Matrix, Ray3D, or sequence')
            else:
                _normal = normal
        else:
            _normal = Matrix(plane.normal_vector)
        mag_incident = sqrt(sum((i ** 2 for i in _incident)))
        mag_normal = sqrt(sum((i ** 2 for i in _normal)))
        mag_refracted = sqrt(sum((i ** 2 for i in refracted)))
        _incident /= mag_incident
        _normal /= mag_normal
        refracted /= mag_refracted
        i = acos(_incident.dot(_normal))
        r = acos(refracted.dot(_normal))
        return i - r

def brewster_angle(medium1, medium2):
    n1 = refractive_index_of_medium(medium1)
    n2 = refractive_index_of_medium(medium2)
    return atan2(n2, n1)

def critical_angle(medium1, medium2):
    n1 = refractive_index_of_medium(medium1)
    n2 = refractive_index_of_medium(medium2)
    if n2 > n1:
        raise ValueError('Total internal reflection impossible for n1 < n2')
    else:
        return asin(n2 / n1)

def lens_makers_formula(n_lens, n_surr, r1, r2, d=0):
    if isinstance(n_lens, Medium):
        n_lens = n_lens.refractive_index
    else:
        n_lens = sympify(n_lens)
    if isinstance(n_surr, Medium):
        n_surr = n_surr.refractive_index
    else:
        n_surr = sympify(n_surr)
    d = sympify(d)
    focal_length = 1 / ((n_lens - n_surr) / n_surr * (1 / r1 - 1 / r2 + (n_lens - n_surr) * d / (n_lens * r1 * r2)))
    if focal_length == zoo:
        return S.Infinity
    return focal_length

def mirror_formula(focal_length=None, u=None, v=None):
    if focal_length and u and v:
        raise ValueError('Please provide only two parameters')
    focal_length = sympify(focal_length)
    u = sympify(u)
    v = sympify(v)
    if u is oo:
        _u = Symbol('u')
    if v is oo:
        _v = Symbol('v')
    if focal_length is oo:
        _f = Symbol('f')
    if focal_length is None:
        if u is oo and v is oo:
            return Limit(Limit(_v * _u / (_v + _u), _u, oo), _v, oo).doit()
        if u is oo:
            return Limit(v * _u / (v + _u), _u, oo).doit()
        if v is oo:
            return Limit(_v * u / (_v + u), _v, oo).doit()
        return v * u / (v + u)
    if u is None:
        if v is oo and focal_length is oo:
            return Limit(Limit(_v * _f / (_v - _f), _v, oo), _f, oo).doit()
        if v is oo:
            return Limit(_v * focal_length / (_v - focal_length), _v, oo).doit()
        if focal_length is oo:
            return Limit(v * _f / (v - _f), _f, oo).doit()
        return v * focal_length / (v - focal_length)
    if v is None:
        if u is oo and focal_length is oo:
            return Limit(Limit(_u * _f / (_u - _f), _u, oo), _f, oo).doit()
        if u is oo:
            return Limit(_u * focal_length / (_u - focal_length), _u, oo).doit()
        if focal_length is oo:
            return Limit(u * _f / (u - _f), _f, oo).doit()
        return u * focal_length / (u - focal_length)

def lens_formula(focal_length=None, u=None, v=None):
    if focal_length and u and v:
        raise ValueError('Please provide only two parameters')
    focal_length = sympify(focal_length)
    u = sympify(u)
    v = sympify(v)
    if u is oo:
        _u = Symbol('u')
    if v is oo:
        _v = Symbol('v')
    if focal_length is oo:
        _f = Symbol('f')
    if focal_length is None:
        if u is oo and v is oo:
            return Limit(Limit(_v * _u / (_u - _v), _u, oo), _v, oo).doit()
        if u is oo:
            return Limit(v * _u / (_u - v), _u, oo).doit()
        if v is oo:
            return Limit(_v * u / (u - _v), _v, oo).doit()
        return v * u / (u - v)
    if u is None:
        if v is oo and focal_length is oo:
            return Limit(Limit(_v * _f / (_f - _v), _v, oo), _f, oo).doit()
        if v is oo:
            return Limit(_v * focal_length / (focal_length - _v), _v, oo).doit()
        if focal_length is oo:
            return Limit(v * _f / (_f - v), _f, oo).doit()
        return v * focal_length / (focal_length - v)
    if v is None:
        if u is oo and focal_length is oo:
            return Limit(Limit(_u * _f / (_u + _f), _u, oo), _f, oo).doit()
        if u is oo:
            return Limit(_u * focal_length / (_u + focal_length), _u, oo).doit()
        if focal_length is oo:
            return Limit(u * _f / (u + _f), _f, oo).doit()
        return u * focal_length / (u + focal_length)

def hyperfocal_distance(f, N, c):
    f = sympify(f)
    N = sympify(N)
    c = sympify(c)
    return 1 / (N * c) * f ** 2

def transverse_magnification(si, so):
    si = sympify(si)
    so = sympify(so)
    return -(si / so)
