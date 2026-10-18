"""Lazy re-exports for compatible import paths over existing implementations.

Resolve at attribute access so packages imported during symox initialization do
not capture unfinished constructors or bypass later native dispatch adapters.
"""
from importlib import import_module


def export(namespace, source, names):
    names = names.split() if isinstance(names, str) else names
    sources = namespace.setdefault('_api_export_sources', {})
    sources.update((name, source) for name in names)
    existing = namespace.get('__all__', [name for name in namespace if not name.startswith('_')])
    namespace['__all__'] = list(dict.fromkeys(existing + list(names)))

    def resolve(name):
        if name in sources:
            return getattr(import_module(sources[name]), name)
        raise AttributeError(f"module {namespace['__name__']!r} has no attribute {name!r}")

    def directory():
        return sorted(set(namespace) | set(sources))

    namespace['__getattr__'] = resolve
    namespace['__dir__'] = directory
