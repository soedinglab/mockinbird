"""
Classes for intervals on genomic sequences
"""


class GenomicInterval:
    """
    Base class for an interval on a genomic sequence

    Parameters
    ----------
    start : int
        the start of the interval
    end : int
        the end of the interval
    one_based : bool
        one_based interval bounds if true else zero-based bounds

    Attributes
    ----------
    start : int
        the start of the interval
    end : int
        the end of the interval
   """

    def __init__(self, start, end, one_based=True):
        self._start = start
        self._end = end
        if not one_based:
            int_start, int_end = start + 1, end
        else:
            int_start, int_end = start, end

        if int_end < int_start:
            raise ValueError("The length of an interval must be positive")
        self._int_start = int_start
        self._int_end = int_end

    def overlaps(self, interval):
        """
        Return `True` if the given interval overlaps this interval
        else `False`

        Parameters
        ----------
        interval : GenomicInterval
            The interval which is checked for an overlap
        """
        return (
            self._int_start <= interval._int_end and
            interval._int_start <= self._int_end
        )

    def __repr__(self):
        return "(%s, %s)" % (self._start, self._end)

    @property
    def start(self):
        """Return the start of the interval"""
        return self._start

    @property
    def end(self):
        """Return the end of the interval"""
        return self._end
