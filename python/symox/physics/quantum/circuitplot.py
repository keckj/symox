"""Matplotlib diagrams of ordered quantum circuits, including measurement symbols."""
from .gate import Gate,OneQubitGate,CGate,CNotGate,SwapGate,XGate
from .._algebra import OperatorExpr


def labeller(n,symbol='q'):return [f'{symbol}_{i}' for i in reversed(range(int(n)))]
class Mz(OneQubitGate):measurement=True;gate_name='Mz';gate_name_latex='M_z'
class Mx(OneQubitGate):measurement=True;gate_name='Mx';gate_name_latex='M_x'

def CreateOneQubitGate(name,latexname=None):return type(str(name)+'Gate',(OneQubitGate,),{'gate_name':str(name),'gate_name_latex':latexname or str(name)})
def CreateCGate(name,latexname=None):
    gate=CreateOneQubitGate(name,latexname)
    def controlled(controls,target):return CGate(tuple(controls),gate(target))
    return controlled


class CircuitPlot:
    def __init__(self,c,nqubits,**kwargs):
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle,Circle
        self.circuit=c;self.nqubits=int(nqubits)
        if self.nqubits<1:raise ValueError('circuit requires at least one qubit')
        if isinstance(c,OperatorExpr):
            if len(c._terms)!=1:raise ValueError('a circuit diagram requires one ordered gate product')
            word,coefficient=next(iter(c._terms.items()))
        elif c==1:word=()
        else:raise TypeError('circuit must be a gate or an ordered gate product')
        gates=list(reversed(word))
        if any(not isinstance(g,Gate) for g in gates):raise TypeError('circuit diagrams accept gate and measurement factors')
        self.ngates=len(gates);self.figure=kwargs.get('figure') or plt.figure(figsize=kwargs.get('figsize',(max(4,len(gates)+1),max(2,self.nqubits*.6))))
        self.axes=kwargs.get('axes') or self.figure.add_subplot(111)
        axes=self.axes;labels=kwargs.get('labels',labeller(self.nqubits));self.gates=gates
        for bit in range(self.nqubits):axes.plot([0,len(gates)+1],[bit,bit],color='black',lw=1)
        for column,gate in enumerate(gates,1):
            targets=gate.targets;controls=getattr(gate,'controls',())
            if any(t>=self.nqubits for t in tuple(targets)+tuple(controls)):raise ValueError('gate target lies outside the diagram')
            if controls:
                if isinstance(gate,CGate):targets=gate.gate.targets
                else:targets=tuple(t for t in targets if t not in controls)
                axes.plot([column,column],[min(tuple(targets)+tuple(controls)),max(tuple(targets)+tuple(controls))],color='black')
                for bit in controls:axes.add_patch(Circle((column,bit),.065,color='black'))
            for bit in targets:
                if isinstance(gate,SwapGate):axes.plot([column-.1,column+.1],[bit-.1,bit+.1],color='black');axes.plot([column-.1,column+.1],[bit+.1,bit-.1],color='black');continue
                if isinstance(gate,CNotGate) or isinstance(gate,CGate) and isinstance(gate.gate,XGate):
                    axes.add_patch(Circle((column,bit),.15,fill=False,color='black'));axes.plot([column,column],[bit-.15,bit+.15],color='black');continue
                label=getattr(gate.gate if isinstance(gate,CGate) else gate,'gate_name',type(gate).__name__.replace('Gate',''))
                axes.add_patch(Rectangle((column-.22,bit-.22),.44,.44,facecolor='white',edgecolor='black',zorder=3));axes.text(column,bit,label,ha='center',va='center',zorder=4)
        axes.set_yticks(list(reversed(range(self.nqubits)))),axes.set_yticklabels(labels)
        axes.set_xlim(-.2,len(gates)+1.2);axes.set_ylim(-.6,self.nqubits-.4);axes.set_xticks([])
        for spine in axes.spines.values():spine.set_visible(False)
    def savefig(self,*args,**kwargs):return self.figure.savefig(*args,**kwargs)

def circuit_plot(c,nqubits,**kwargs):return CircuitPlot(c,nqubits,**kwargs)
__all__=['CircuitPlot','circuit_plot','labeller','Mz','Mx','CreateOneQubitGate','CreateCGate']
