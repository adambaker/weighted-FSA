from fsa.nat_class_set import NaturalClassSet

class TestNatClassSet:
    def setup_method(self, method):
        self.alphabet = set('aeuptkbdszh')
        self.classes = set(
            [
                frozenset('aeu'),
                frozenset('ptkbdzsh'),
                frozenset('ptkbd'),
                frozenset('szh'),
                frozenset('bdz'),
                frozenset('ptksh'),
                frozenset('ptbdsz')
            ])
        self.nat_classes = NaturalClassSet(self.alphabet, self.classes)

    def classes_in_nat_class_set(self, classes):
        for class_ in classes:
            assert frozenset(class_) in self.nat_classes
    
    def test_class_contains_singletons(self):
        self.classes_in_nat_class_set(self.alphabet)

    def test_class_has_natural_classes(self):
        self.classes_in_nat_class_set(self.classes)

    def test_intersection_classes_in_class(self):
        self.classes_in_nat_class_set(
            ['bd', 'sh', 'pts', 'pt', 'ptdb', 'ptk', 'sz'])
      
    def test_class_size(self):
        assert len(self.nat_classes) == 25

from fsa.arc_set import NatClassArcSet

class TestNatClassArcSet:
    def setup_method(self, method):
        self.coronal = frozenset('tpbdrzlnm')
        self.stop = frozenset('bd')
        self.nasal = frozenset('nm')
        self.classes = NaturalClassSet(
            'ptkbdgrszlmn', set([self.coronal, self.nasal, self.stop]))
        
    def test_simple_arc_set(self):
        arcset = NatClassArcSet('1', '1', 'tpbdrzlnm', 3, self.classes)
        assert arcset.labels == set([self.coronal])

    def test_complex_arc_set(self):
        arcset = NatClassArcSet('1', '1', 'bdrzlmn', 2, self.classes)
        labels = arcset.labels
        
        assert len(labels) == 5
        assert self.stop in labels
        assert self.nasal in labels
        assert frozenset('z') in labels
        assert frozenset('l') in labels
        assert frozenset('r') in labels
