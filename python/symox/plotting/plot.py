"""Sample native expressions and render through an optional Matplotlib backend.

Sampling is deterministic on a uniform grid. No external symbolic engine is
used. Matplotlib is imported only when a figure is shown or saved.
"""
import math
import operator
import numpy as np


def _s():
    import symox
    return symox


def _range(value):
    return isinstance(value,(tuple,list)) and len(value)==3 and getattr(value[0],'is_Symbol',False)


def _count(value):
    try:
        n=operator.index(value)
    except TypeError:
        raise ValueError('Sample counts must be integers') from None
    if not 2 <= n <= 100000:
        raise ValueError('Sample counts must be between 2 and 100000')
    return n


def _grid(bound,n):
    symbol,start,end=bound
    start,end=float(start),float(end)
    if not math.isfinite(start) or not math.isfinite(end) or start==end:
        raise ValueError('Plot bounds must be distinct finite real numbers')
    return np.linspace(start,end,_count(n))


def _evaluate(expr,symbols,values):
    function=expr if callable(expr) and not isinstance(expr,_s()._core.Expr) else _s().lambdify(symbols,expr,modules='numpy')
    with np.errstate(all='ignore'):
        result=np.asarray(function(*values),dtype=complex)
    result=np.broadcast_to(result,np.broadcast_shapes(*(np.shape(v) for v in values)))
    real=result.real.copy()
    real[(np.abs(result.imag)>1e-10) | ~np.isfinite(result)]=np.nan
    return real


class BaseSeries:
    """An expression, its parameter ranges and sampled numeric data."""
    is_2Dline=False
    is_3Dline=False
    is_3Dsurface=False
    is_contour=False
    is_implicit=False
    def __init__(self,expressions,ranges,kind,label=None,**options):
        self.expr=expressions[0] if len(expressions)==1 else tuple(expressions)
        self.expressions=tuple(expressions)
        self.ranges=tuple(tuple(r) for r in ranges)
        self.kind=kind
        self.label=str(self.expr) if label is None else str(label)
        self.line_color=options.get('line_color')
        self.surface_color=options.get('surface_color')
        self.rendering_kw=dict(options.get('rendering_kw',{}))
        self.is_2Dline=kind in ('line','parametric')
        self.is_3Dline=kind=='parametric3d'
        self.is_3Dsurface=kind in ('surface','surface_parametric')
        self.is_contour=kind=='contour'
        self.is_implicit=kind=='implicit'
        self.is_parametric=kind in ('parametric','parametric3d','surface_parametric')
        self.adaptive=False
        self.sampling='uniform'
        default=100 if len(ranges)==2 else 1000
        n=options.get('n',default)
        if isinstance(n,(tuple,list)):
            if len(n)!=len(ranges):raise ValueError('One sample count is required per range')
            self.n=tuple(_count(v) for v in n)
        else:
            self.n=(_count(n),)*len(ranges)
        self.n=tuple(_count(options.get('n%d'%(i+1),v)) for i,v in enumerate(self.n))
        self.detect_poles=bool(options.get('detect_poles',False))
        self._implicit_mode=options.get('_implicit_mode','contour')
        self._implicit_relation=options.get('_implicit_relation')
        self._params=options.get('params',{})

    def get_label(self,use_latex=False):
        return self.label

    def get_data(self):
        symbols=[r[0] for r in self.ranges]
        grids=[_grid(r,n) for r,n in zip(self.ranges,self.n)]
        values=grids if len(grids)==1 else np.meshgrid(*grids,indexing='xy')
        expressions=[e.subs(self._params) if self._params and hasattr(e,'subs') else e for e in self.expressions]
        sampled=[_evaluate(e,symbols,values) for e in expressions]
        if self.kind=='line':
            y=sampled[0]
            if self.detect_poles and len(y)>2:
                # Break opposite-sign jumps that exceed the local vertical scale.
                delta=np.abs(np.diff(y))
                finite=delta[np.isfinite(delta)]
                scale=np.median(finite) if finite.size else 0
                jumps=(np.signbit(y[:-1])!=np.signbit(y[1:])) & (delta>max(1,100*scale))
                y[np.flatnonzero(jumps)+1]=np.nan
            return grids[0],y
        if self.kind in ('parametric','parametric3d','surface_parametric'):
            return (*sampled,*values) if self.is_parametric and self.kind!='surface_parametric' else tuple(sampled)
        if self.kind=='implicit':
            z=sampled[0]
            if self._implicit_relation:
                name=self._implicit_relation
                z={'LessThan':lambda:z<=0,'StrictLessThan':lambda:z<0,
                   'GreaterThan':lambda:z>=0,'StrictGreaterThan':lambda:z>0,
                   'Unequality':lambda:z!=0}[name]().astype(float)
            return *values,z,self._implicit_mode
        return *values,sampled[0]


class MatplotlibBackend:
    """Lazily create a Matplotlib figure for a Plot or PlotGrid."""
    def __init__(self,parent):
        self.parent=parent
        self.fig=None
        self.ax=None

    def _draw(self,ax,series):
        data=series.get_data()
        kw=dict(series.rendering_kw)
        if series.is_2Dline or series.is_3Dline:
            if series.line_color is not None:kw.setdefault('color',series.line_color)
            kw.setdefault('label',series.label)
            ax.plot(*data[:3 if series.is_3Dline else 2],**kw)
        elif series.is_3Dsurface:
            if series.surface_color is not None:kw.setdefault('color',series.surface_color)
            ax.plot_surface(*data[:3],**kw)
        elif series.is_contour:
            ax.contour(*data[:3],**kw)
        elif series.is_implicit:
            x,y,z,mode=data
            kw.setdefault('colors',series.line_color or 'blue')
            if mode=='contour':ax.contour(x,y,z,levels=[0],**kw)
            else:ax.contourf(x,y,z,levels=[0.5,1.5],**kw)

    def _axes(self,ax,plot):
        for series in plot:
            self._draw(ax,series)
        options=plot.options
        for key in ('title','xlabel','ylabel','zlabel','xlim','ylim','zlim','xscale','yscale'):
            if key in options and hasattr(ax,'set_'+key):getattr(ax,'set_'+key)(options[key])
        if options.get('legend'):ax.legend()
        if 'axis' in options and not options['axis']:ax.set_axis_off()

    def process_series(self):
        if self.fig is not None:return
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError('Plot rendering requires matplotlib; install symox[plotting]') from None
        parent=self.parent
        self.fig=plt.figure(figsize=parent.size)
        try:
            if isinstance(parent,PlotGrid):
                axes=[]
                for i,plot in enumerate(parent._plots):
                    projection='3d' if any(s.is_3Dline or s.is_3Dsurface for s in plot) else None
                    ax=self.fig.add_subplot(parent.nrows,parent.ncolumns,i+1,projection=projection)
                    self._axes(ax,plot)
                    axes.append(ax)
                self.ax=axes
            else:
                projection='3d' if any(s.is_3Dline or s.is_3Dsurface for s in parent) else None
                self.ax=self.fig.add_subplot(111,projection=projection)
                self._axes(self.ax,parent)
            self.fig.tight_layout()
        except Exception:
            self.close()
            raise

    def show(self):
        self.process_series()
        import matplotlib.pyplot as plt
        plt.show()

    def save(self,path,**kwargs):
        self.process_series()
        self.fig.savefig(path,**kwargs)

    def close(self):
        if self.fig is not None:
            import matplotlib.pyplot as plt
            plt.close(self.fig)
            self.fig=None
            self.ax=None


class TextBackend:
    def __init__(self,parent):self.parent=parent
    def show(self):
        for series in self.parent:
            if series.kind!='line':raise NotImplementedError('The text backend supports scalar 2D lines')
            textplot(series.expr,*series.ranges[0][1:])
    def save(self,path,**kwargs):
        from contextlib import redirect_stdout
        with open(path,'w') as stream,redirect_stdout(stream):self.show()
    def close(self):pass


plot_backends={'matplotlib':MatplotlibBackend,'default':MatplotlibBackend,'text':TextBackend}


class Plot:
    def __init__(self,*series,**options):
        self._series=list(series)
        self.options=dict(options)
        self.size=options.get('size')
        backend=options.get('backend','default')
        self.backend=plot_backends[backend] if isinstance(backend,str) else backend
        self._backend=self.backend(self)
    def __len__(self):return len(self._series)
    def __iter__(self):return iter(self._series)
    def __getitem__(self,index):return self._series[index]
    def __setitem__(self,index,value):
        if not isinstance(value,BaseSeries):raise TypeError('Expected a plot series')
        self._series[index]=value
        self._backend.close()
    def __delitem__(self,index):
        del self._series[index]
        self._backend.close()
    def append(self,value):
        if not isinstance(value,BaseSeries):raise TypeError('Expected a plot series')
        self._series.append(value)
        self._backend.close()
    def extend(self,value):
        for series in value:self.append(series)
    def show(self):return self._backend.show()
    def save(self,path,**kwargs):return self._backend.save(path,**kwargs)
    def close(self):return self._backend.close()
    @property
    def fig(self):
        self._backend.process_series()
        return self._backend.fig
    @property
    def ax(self):
        self._backend.process_series()
        return self._backend.ax


class PlotGrid:
    def __init__(self,nrows,ncolumns,*args,show=True,size=None,**kwargs):
        self.nrows,self.ncolumns=operator.index(nrows),operator.index(ncolumns)
        if self.nrows<1 or self.ncolumns<1 or len(args)>self.nrows*self.ncolumns:
            raise ValueError('Plot grid dimensions must be positive and fit all plots')
        if not all(isinstance(p,Plot) for p in args):raise TypeError('PlotGrid requires Plot objects')
        self._plots=list(args)
        self.size=size
        self._backend=MatplotlibBackend(self)
        if show:self.show()
    def __getitem__(self,index):return self._plots[index]
    @property
    def fig(self):
        self._backend.process_series()
        return self._backend.fig
    def show(self):return self._backend.show()
    def save(self,path,**kwargs):return self._backend.save(path,**kwargs)
    def close(self):return self._backend.close()


def _groups(args,arity,dimensions):
    s=_s()
    if not args:raise ValueError('At least one expression is required')
    global_ranges=[a for a in args if _range(a)]
    entries=[a for a in args if not _range(a)]
    if len(global_ranges)>dimensions:raise ValueError('Too many plotting ranges')
    if arity>1 and len(entries)==arity and not any(isinstance(a,(tuple,list)) for a in entries):
        entries=[entries]
    for entry in entries:
        label=None
        local=[]
        if isinstance(entry,(tuple,list)):
            exprs=[]
            for part in entry:
                if _range(part):local.append(part)
                elif isinstance(part,str) and len(exprs)>=arity:label=part
                else:exprs.append(part)
        else:exprs=[entry]
        if len(exprs)!=arity:raise ValueError('Expected %d expression(s) per series'%arity)
        exprs=[e if callable(e) else s.sympify(e) for e in exprs]
        ranges=local or list(global_ranges)
        if len(ranges)>dimensions:raise ValueError('Too many plotting ranges')
        symbols=set().union(*(getattr(e,'free_symbols',set()) for e in exprs))
        for r in ranges:symbols.discard(r[0])
        for sym in sorted(symbols,key=str):
            if len(ranges)==dimensions:raise ValueError('Unbound symbols in plotting expression')
            ranges.append((sym,-10,10))
        while len(ranges)<dimensions:ranges.append((s.Dummy('t%d'%len(ranges)),-10,10))
        if len({r[0] for r in ranges})!=dimensions:raise ValueError('Plot range symbols must be distinct')
        for r in ranges:_grid(r,2)
        yield exprs,ranges,label


def _plot(args,kind,arity,dimensions,show,options):
    params=options.get('params',{})
    if params:
        args=tuple(a.subs(params) if hasattr(a,'subs') and not _range(a) else a for a in args)
    series=[BaseSeries(e,r,kind,label,**options) for e,r,label in _groups(args,arity,dimensions)]
    result=Plot(*series,**options)
    if show:result.show()
    return result


def plot(*args,show=True,**kwargs):
    return _plot(args,'line',1,1,show,kwargs)


def plot_parametric(*args,show=True,**kwargs):
    return _plot(args,'parametric',2,1,show,kwargs)


def plot3d(*args,show=True,**kwargs):
    return _plot(args,'surface',1,2,show,kwargs)


def plot3d_parametric_line(*args,show=True,**kwargs):
    return _plot(args,'parametric3d',3,1,show,kwargs)


def plot3d_parametric_surface(*args,show=True,**kwargs):
    return _plot(args,'surface_parametric',3,2,show,kwargs)


def plot_contour(*args,show=True,**kwargs):
    return _plot(args,'contour',1,2,show,kwargs)


def plot_implicit(expr,x_var=None,y_var=None,adaptive=True,depth=0,n=300,line_color='blue',show=True,**kwargs):
    s=_s();expr=s.sympify(expr)
    mode='contour'
    relation=None
    if getattr(expr,'is_Relational',False):
        name=expr.class_name
        if name!='Equality':mode='contourf';relation=name
        expr=expr.lhs-expr.rhs
    elif expr.class_name in ('And','Or','Not'):
        raise NotImplementedError('Compound implicit regions are not implemented')
    kwargs.update(n=n,line_color=line_color,_implicit_mode=mode,_implicit_relation=relation)
    ranges=tuple(v for v in (x_var,y_var) if v is not None)
    return _plot((expr,*ranges),'implicit',1,2,show,kwargs)


def textplot(expr,a,b,W=55,H=21):
    """Print an ASCII graph; a,b may be endpoints or (symbol,start,end)."""
    W,H=operator.index(W),operator.index(H)
    if W<2 or H<2:raise ValueError('Text plot dimensions must be at least two')
    s=_s();expr=s.sympify(expr)
    if isinstance(a,_core_expr_type()) and a.is_Symbol:
        raise TypeError('textplot expects numeric endpoints a and b')
    symbols=sorted(expr.free_symbols,key=str)
    if len(symbols)>1:raise ValueError('Text plots require at most one free symbol')
    symbol=symbols[0] if symbols else s.Dummy('x')
    x,y=BaseSeries([expr],[(symbol,a,b)],'line',n=W).get_data()
    valid=np.isfinite(y)
    if not np.any(valid):raise ValueError('No finite real samples in plotting range')
    low,high=np.min(y[valid]),np.max(y[valid])
    if high==low:low-=1;high+=1
    canvas=[[' ']*W for _ in range(H)]
    if low<=0<=high:
        row=round((high/(high-low))*(H-1))
        canvas[row]=['-']*W
    if min(a,b)<=0<=max(a,b):
        col=round((0-float(a))/(float(b)-float(a))*(W-1))
        for row in canvas:row[col]='|'
    for col,value in enumerate(y):
        if np.isfinite(value):canvas[round((high-value)/(high-low)*(H-1))][col]='*'
    print('\n'.join(''.join(row) for row in canvas))


def _core_expr_type():
    return _s()._core.Expr


def PygletPlot(*args,**kwargs):
    """Legacy plot entry point backed by Matplotlib, including 3D surfaces.

    This provides composable figures instead of Pyglet's OpenGL window.
    """
    dimensions=sum(_range(a) for a in args)
    return (plot3d if dimensions==2 else plot)(*args,**kwargs)


def check_arguments(args,expr_len,nb_of_free_symbols):
    if not args:return []
    return [tuple(expressions)+tuple(ranges)+((label,) if label is not None else ()) for expressions,ranges,label in _groups(args,expr_len,nb_of_free_symbols)]


def plot_factory(*args,**kwargs):
    """Construct a backend-compatible plot container from existing series."""
    backend=kwargs.pop('backend','default')
    if isinstance(backend,str):
        if backend not in ('default','matplotlib','text'):raise ValueError('unknown plot backend '+backend)
        return Plot(*args,backend=backend,**kwargs)
    if isinstance(backend,type) and issubclass(backend,Plot):return backend(*args,**kwargs)
    raise TypeError('backend must be a known name or a Plot subclass')
