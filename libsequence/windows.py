from .polytable import SimData
class Windows:
    """
    An iterable list of sliding windows created from a :class:`libsequence.polytable.PolyTable`
    """
    def __init__(self, pt, window_size, step_len, starting_pos = 0.,  ending_pos = 1):
        self.windows=[]
        if isinstance(pt,SimData):
            from .windows_cpp import SimDataWindows
            temp = SimDataWindows(pt,window_size,step_len,starting_pos,ending_pos)
        else:
            from .windows_cpp import PolySitesWindows
            temp = PolySitesWindows(pt,window_size,step_len,starting_pos,ending_pos)
        for i in range(len(temp)):
            self.windows.append(temp[i])
        temp=None
    def __iter__(self):
        return iter(self.windows)
    def __next__(self):
        return next(self.windows)
    def __getitem__(self,i):
        return self.windows[i]
    def __len__(self):
        return len(self.windows)
