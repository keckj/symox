"""Objects representing generated program statements."""

from .ast import Assignment, CodeBlock, Comment

__all__ = ['Assignment', 'CodeBlock', 'Comment']
from .ast import (Attribute,Declaration,For,FunctionCall,FunctionDefinition,FunctionPrototype,Print,Scope,Variable,While,aug_assign)
__all__ += ['Attribute','Declaration','For','FunctionCall','FunctionDefinition','FunctionPrototype','Print','Scope','Variable','While','aug_assign']
