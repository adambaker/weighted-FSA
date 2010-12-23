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
        self.labels = _get_covering_labels(frozenset(letters),
             nat_class_set.classes)


def _get_covering_labels( letters, classes ):
    ''' 
    >>> from nat_class_set import NaturalClassSet
    >>> coronal = frozenset('tpbdrzlnm')
    >>> stop = frozenset('bd')
    >>> nasal = frozenset('nm')
    >>> classes = NaturalClassSet('ptkbdgrszlmn', set([coronal, nasal, stop]))
    >>> labels = _get_covering_labels(set('bdrzlmn'), classes.classes)
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
    >>> labels = _get_covering_labels(set('tpbdrzlnm'), classes.classes)
    >>> labels == set([coronal])
    True
    '''  #this doctest has been adapted to test NatClassArcSet in 
         #test/test_nat_class_set.py.
    coverage = set([])
    covering_labels = set([])
    classes = set(classes) 
    while True:
        to_cover = letters.difference(coverage)
        to_add = frozenset([])
        for class_ in classes:
            if class_.issubset(to_cover) and len(class_) > len(to_add):
                to_add = class_
        if len(to_add) > 0:
            covering_labels.add(to_add)
            coverage.update(to_add)
            classes.remove(to_add)
        if len(to_add) <= 1:break
    #only singletons left to fill out the rest of the letters.
    for letter in letters.difference(coverage):
        covering_labels.add(frozenset([letter]))
    return frozenset(covering_labels)

if __name__ == '__main__':
    import doctest
    doctest.testmod()

