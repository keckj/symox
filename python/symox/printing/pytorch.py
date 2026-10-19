"""SymPy-compatible PyTorch printer import path.

Symox-specific tensor execution and source exporters live in symox.extra.
"""
from ..extra.printing.pytorch import TorchPrinter, torchcode as torch_code
__all__ = ['TorchPrinter', 'torch_code']
