"""Incremental memoization for ordinary and associated recurrence sequences."""
from functools import wraps
from threading import RLock
import operator

def recurrence_memo(initial):
    cache=list(initial);lock=RLock()
    def decorate(function):
        @wraps(function)
        def sequence(n):
            n=operator.index(n)
            with lock:
                # As in the reference, an index below the cache length,
                # including a negative one, reads the cache directly
                # (negative indices count from its end).
                if n<len(cache):return cache[n]
                for i in range(len(cache),n+1):cache.append(function(i,cache))
                return cache[n]
        sequence.cache_length=lambda:len(cache)
        sequence.fetch_item=lambda index:cache[index]
        return sequence
    return decorate

def assoc_recurrence_memo(base_seq):
    cache=[];lock=RLock()
    def decorate(function):
        @wraps(function)
        def sequence(n,m):
            n,m=operator.index(n),operator.index(m)
            with lock:
                # The reference indexes the cache rows directly: negative
                # n and m count from the end, m > n raises IndexError.
                if n<len(cache):return cache[n][m]
                for i in range(len(cache),n+1):
                    cache.append([base_seq(i)])
                    for j in range(1,i+1):cache[i].append(function(i,j,cache))
                return cache[n][m]
        return sequence
    return decorate
__all__=['recurrence_memo','assoc_recurrence_memo']
