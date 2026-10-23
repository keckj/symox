"""Deterministic finite diagram layouts and executable Xy-pic output."""
from collections import defaultdict,deque
from .baseclasses import Diagram,Object,IdentityMorphism,CompositeMorphism


class DiagramGrid:
    def __init__(self,diagram,groups=None,**hints):
        if not isinstance(diagram,Diagram):
            raise TypeError('expected a Diagram')
        unknown = set(hints)-{'layout','transpose'}
        if unknown:
            raise NotImplementedError('unsupported layout options: '+', '.join(sorted(unknown)))
        if hints.get('layout','generic') not in ('generic','sequential'):
            raise ValueError('layout must be generic or sequential')
        self.diagram = diagram
        merged = dict(diagram.premises)
        for m,props in diagram.conclusions.items():
            merged[m] = merged.get(m,frozenset()) | props
        self._morphisms = {m:p for m,p in merged.items() if not isinstance(m,IdentityMorphism)
                           and (not isinstance(m,CompositeMorphism) or p)}
        objects = sorted(diagram.objects,key=lambda o:o.name)
        if groups is not None:
            ordered = []
            for group in groups:
                group = (group,) if isinstance(group,Object) else group
                for obj in sorted(group,key=lambda o:o.name):
                    if obj not in objects or obj in ordered:
                        raise ValueError('groups must contain distinct diagram objects')
                    ordered.append(obj)
            objects = ordered+[o for o in objects if o not in ordered]
        depth = {o:0 for o in objects}
        incoming = {o:0 for o in objects}
        outgoing = defaultdict(set)
        for m in self._morphisms:
            if m.domain != m.codomain and m.codomain not in outgoing[m.domain]:
                outgoing[m.domain].add(m.codomain)
                incoming[m.codomain] += 1
        queue = deque(o for o in objects if incoming[o] == 0)
        processed = set()
        while queue:
            obj = queue.popleft()
            processed.add(obj)
            for target in sorted(outgoing[obj],key=lambda o:o.name):
                depth[target] = max(depth[target],depth[obj]+1)
                incoming[target] -= 1
                if incoming[target] == 0:
                    queue.append(target)
        # Cyclic components get distinct cells after the acyclic layout.
        column = max((depth[o] for o in processed),default=-1)+1
        for obj in objects:
            if obj not in processed:
                depth[obj] = column
                column += 1
        layers = defaultdict(list)
        for obj in objects:
            layers[depth[obj]].append(obj)
        width = max(depth.values(),default=-1)+1
        height = max((len(v) for v in layers.values()),default=0)
        self._grid = [[None]*width for _ in range(height)]
        for col,layer in layers.items():
            for row,obj in enumerate(layer):
                self._grid[row][col] = obj
        if hints.get('transpose',False):
            self._grid = [list(row) for row in zip(*self._grid)]

    width = property(lambda self:len(self._grid[0]) if self._grid else 0)
    height = property(lambda self:len(self._grid))
    morphisms = property(lambda self:dict(self._morphisms))

    def __getitem__(self,indices):
        i,j = indices
        return self._grid[i][j]

    def __str__(self):
        return '\n'.join(' '.join(o.name if o is not None else '.' for o in row) for row in self._grid)


class ArrowStringDescription:
    def __init__(self,unit='mm',curving='',curving_amount=0,looping_start='',looping_end='',
                 horizontal_direction='',vertical_direction='',label_position='^',label=''):
        self.unit,self.curving,self.curving_amount = unit,curving,curving_amount
        self.looping_start,self.looping_end = looping_start,looping_end
        self.horizontal_direction,self.vertical_direction = horizontal_direction,vertical_direction
        self.label_position,self.label = label_position,label
        self.label_displacement,self.arrow_style = '',''
        self.forced_label_position = False

    def __str__(self):
        curve = '@/%s%s%s/' % (self.curving,self.curving_amount,self.unit) if self.curving else ''
        loop = '@(%s,%s)' % (self.looping_start,self.looping_end) if self.looping_start and self.looping_end else ''
        style = '@'+self.arrow_style if self.arrow_style else ''
        return r'\ar%s%s%s[%s%s]%s%s{%s}' % (curve,loop,style,self.horizontal_direction,self.vertical_direction,self.label_position,self.label_displacement,self.label)


def _morphism_label(morphism,latex):
    from .. import Symbol
    from .baseclasses import CompositeMorphism,IdentityMorphism,NamedMorphism
    if isinstance(morphism,IdentityMorphism):
        return r'\mathrm{id}_{%s}' % latex(morphism.domain)
    if isinstance(morphism,CompositeMorphism):
        return r'\circ '.join(latex(Symbol(c.name)) for c in reversed(morphism.components))
    if isinstance(morphism,NamedMorphism):
        return latex(Symbol(morphism.name))
    return ''


class XypicDiagramDrawer:
    def __init__(self):
        self.unit = 'mm'
        self.default_curving_amount = 3
        self.default_curving_step = 4
        self.arrow_formatters = {}
        self.default_arrow_formatter = None

    def draw(self,diagram,grid,masked=None,diagram_format=''):
        from ..printing import latex
        if not isinstance(diagram,Diagram) or not isinstance(grid,DiagramGrid):
            raise TypeError('expected a Diagram and DiagramGrid')
        positions = {grid[i,j]:(i,j) for i in range(grid.height) for j in range(grid.width) if grid[i,j] is not None}
        if set(positions) != set(diagram.objects):
            raise ValueError('grid does not contain exactly the diagram objects')
        arrows = defaultdict(list)
        parallel = defaultdict(int)
        for m,props in sorted(grid.morphisms.items(),key=lambda item:repr(item[0])):
            if masked is not None and m in masked:
                continue
            i,j = positions[m.domain]
            ti,tj = positions[m.codomain]
            desc = ArrowStringDescription(unit=self.unit,label=_morphism_label(m,latex))
            desc.horizontal_direction = ('r' if tj >= j else 'l')*abs(tj-j)
            desc.vertical_direction = ('d' if ti >= i else 'u')*abs(ti-i)
            key = frozenset((m.domain,m.codomain))
            count = parallel[key]
            parallel[key] += 1
            if m.domain == m.codomain:
                desc.looping_start,desc.looping_end = 'r','u'
            elif count:
                desc.curving = '^' if count % 2 else '_'
                desc.curving_amount = self.default_curving_amount+((count-1)//2)*self.default_curving_step
            if self.default_arrow_formatter is not None:
                self.default_arrow_formatter(desc)
            for prop in sorted(props,key=str):
                formatter = self.arrow_formatters.get(str(prop))
                if formatter is not None:
                    formatter(desc)
            arrows[m.domain].append(str(desc))
        rows = []
        for i in range(grid.height):
            cells = []
            for j in range(grid.width):
                obj = grid[i,j]
                cells.append('' if obj is None else latex(obj)+' '+ ' '.join(arrows[obj]))
            rows.append(' & '.join(cells).rstrip())
        return '\\xymatrix'+diagram_format+'{\n'+(' \\\\\n'.join(rows))+'\n}'


def xypic_draw_diagram(diagram,masked=None,diagram_format='',groups=None,**hints):
    return XypicDiagramDrawer().draw(diagram,DiagramGrid(diagram,groups,**hints),masked,diagram_format)


def preview_diagram(diagram,masked=None,diagram_format='',groups=None,output='png',viewer=None,euler=True,**hints):
    from ..printing import preview
    preview_options = {key:hints.pop(key) for key in ('filename','outputbuffer','dvioptions') if key in hints}
    source = xypic_draw_diagram(diagram,masked,diagram_format,groups,**hints)
    return preview('$'+source+'$',output=output,viewer=viewer,euler=euler,
                   packages=('xy',),extra_preamble='\n\\xyoption{all}\n',**preview_options)
