class SequenceSample:
    def __init__(self,name,seq):
        self.name = name
        self.seq = seq

    def __eq__(self,other):
        return self.seq == other.seq and self.name == other.name

    def __hash__(self):
        return hash((self.name,self.seq))

    def copy(self):
        return type(self)(self.name,self.seq)

    def deepcopy(self):
        return type(self)(self.name,self.seq)
