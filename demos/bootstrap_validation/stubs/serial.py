# stub: PMOIRED imports `serial` at module level but only its optional
# microprogress (Arduino progress bar) helper actually uses it.
class Serial:
    def __init__(self, *a, **k): raise RuntimeError("serial stub")
def tools(*a, **k): raise RuntimeError("serial stub")
