# A collection of arcs with the same source and destination states and weights
# but different letters.

class ArcSet(object):
    def __init__(self, source, dest, letters, weight):
        self.source = source
        self.dest = dest
        self.letters = frozenset(letters)
        self.weight = weight

    def __contains__(self, letter):
        return letter in self.letters

    def __iter__(self):
        return self.letters.__iter__()

    def __len__(self):
        return len(self.letters)

class NatClassArcSet(ArcSet):
    def __init__(self, source, dest, letters, weight, nat_class_set):
        self.source = source
        self.dest = dest
        self.letters = frozenset(letters)
        self.weight = weight
        self.labels = get_covering_labels(frozenset(letters),
             nat_class_set.classes)


def get_covering_labels( letters, classes ):
    ''' 
    >>> from nat_class_set import NaturalClassSet
    >>> coronal = frozenset('tpbdrzlnm')
    >>> stop = frozenset('bd')
    >>> nasal = frozenset('nm')
    >>> classes = NaturalClassSet('ptkbdgrszlmn', set([coronal, nasal, stop]))
    >>> labels = get_covering_labels(set('bdrzlmn'), classes.classes)
    >>> len(labels)
    5
    >>> stop in labels
    True
    >>> nasal in labels
    True
    >>> frozenset('z') in labels
    True
    >>> frozenset('l') in labels
    True
    >>> frozenset('r') in labels
    True
    >>> labels = get_covering_labels(set('tpbdrzlnm'), classes.classes)
    >>> labels == set([coronal])
    True
    '''
    coverage = set([])
    covering_labels = set([])
    classes = classes.union(covering_labels) #making a copy of classes
    while True:
        to_cover = letters.difference(coverage)
        to_add = frozenset([])
        for cls in classes:
            if cls.issubset(to_cover) and \
                            len(cls) > len(to_add):
                to_add = cls
        if len(to_add) > 0:
            covering_labels.add(to_add)
            coverage.update(to_add)
            classes.remove(to_add)
        if len(to_add) <= 1:break
    #only singletons left to fill out the rest of the letters.
    for letter in letters.difference(coverage):
        covering_labels.add(frozenset([letter]))
    return covering_labels

if __name__ == '__main__':
    import doctest
    doctest.testmod()

