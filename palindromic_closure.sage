#*****************************************************************************
#       Copyright (C) 2011 Alexandre Blondin Masse <ablondin@uqac.ca>,
#
#  Distributed under the terms of the GNU General Public License
#  version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

def _iterated_right_palindromic_closure_iterator(self, f=None):
    r"""
    Returns an iterator over the iterated (`f`-)palindromic closure of self.
                        
    INPUT:

    -  ``f`` - involution (default: None) on the alphabet of self. It must 
       be callable on letters as well as words (e.g. WordMorphism).
                         
    OUTPUT:

        iterator -- the iterated (`f`-)palindromic closure of self
        
    EXAMPLES::

        sage: w = Word('abc')
        sage: it = w._iterated_right_palindromic_closure_iterator()
        sage: Word(it)
        word: abacaba

    ::

        sage: w = Word('aaa')
        sage: it = w._iterated_right_palindromic_closure_iterator()
        sage: Word(it)
        word: aaa

    ::

        sage: w = Word('abbab')
        sage: it = w._iterated_right_palindromic_closure_iterator()
        sage: Word(it)
        word: ababaabababaababa

    An infinite word::

        sage: t = words.ThueMorseWord('ab')
        sage: it = t._iterated_right_palindromic_closure_iterator()
        sage: Word(it)
        word: ababaabababaababaabababaababaabababaabab...

    TESTS:

    The empty word::

        sage: w = Word()
        sage: it = w._iterated_right_palindromic_closure_iterator()
        sage: it.next()
        Traceback (most recent call last):
        ...
        StopIteration

    REFERENCES:

    -   [1] A. de Luca, A. De Luca, Pseudopalindrome closure operators
        in free monoids, Theoret. Comput. Sci. 362 (2006) 282--300.      
    """
    par = self.parent()
    w = self[:0]
    for letter in self:
        length_before = w.length()
        w = (w*par([letter])).palindromic_closure(f=f)
        length_after = w.length()
        d = length_after - length_before
        for a in w[-d:]:
            yield a

def _iterated_right_palindromic_closure_recursive_iterator(self, f=None):
    r"""
    Returns an iterator over the iterated (`f`-)palindromic closure of
    self.

    INPUT:

    -  ``f`` - involution (default: None) on the alphabet of self.
       It must be callable on letters as well as words
       (e.g. WordMorphism).

    OUTPUT:

        iterator -- the iterated (`f`-)palindromic closure of self

    ALGORITHM:

        For the case of palindromes only, it has been shown in [2] that
        the iterated right palindromic closure of a given word `w`,
        denoted by `IRPC(w)`, may be obtained as follows.
        Let `w` be any word and `x` be a letter. Then

        #. If `x` does not occur in `w`,
           `IRPC(wx) = IRPC(w) \cdot x \cdot IRPC(w)`
        #. Otherwise, write `w = w_1xw_2` such that `x` does not
           occur in `w_2`. Then `IRPC(wx) = IRPC(w) \cdot IRPC(w_1)^{-1}
           \cdot IRPC(w)`

        This formula is directly generalized to the case of
        `f`-palindromes. See [1] for more details.
        
    EXAMPLES::

        sage: w = Word('abc')
        sage: it = w._iterated_right_palindromic_closure_recursive_iterator()
        sage: Word(it)
        word: abacaba

    ::

        sage: w = Word('aaa')
        sage: it = w._iterated_right_palindromic_closure_recursive_iterator()
        sage: Word(it)
        word: aaa

    ::

        sage: w = Word('abbab')
        sage: it = w._iterated_right_palindromic_closure_recursive_iterator()
        sage: Word(it)
        word: ababaabababaababa

    An infinite word::

        sage: t = words.ThueMorseWord('ab')
        sage: it = t._iterated_right_palindromic_closure_recursive_iterator()
        sage: Word(it)
        word: ababaabababaababaabababaababaabababaabab...

    TESTS:

    The empty word::

        sage: w = Word()
        sage: it = w._iterated_right_palindromic_closure_recursive_iterator()
        sage: it.next()
        Traceback (most recent call last):
        ...
        StopIteration

    REFERENCES:

    -   [1] A. de Luca, A. De Luca, Pseudopalindrome closure operators
        in free monoids, Theoret. Comput. Sci. 362 (2006) 282--300.      
    -   [2] J. Justin, Episturmian morphisms and a Galois theorem on
        continued fractions, RAIRO Theoret. Informatics Appl. 39 (2005)
        207-215.
    """
    parent = self.parent()
    ipcw = self[:0]
    lengths = []
    for i, letter in enumerate(self):
        lengths.append(ipcw.length())
        w = self[:i]
        pos = w.rfind(parent([letter]))
        if pos == -1:
            to_append = parent([letter]).palindromic_closure(f=f) + ipcw
        else:
            to_append = ipcw[lengths[pos]:]
        ipcw += to_append
        for a in to_append:
            yield a

def iterated_right_palindromic_closure(self, f=None, algorithm='recursive'):
    r"""
    Returns the iterated (`f`-)palindromic closure of self.
                        
    INPUT:

    -  ``f`` - involution (default: None) on the alphabet of self. It must 
       be callable on letters as well as words (e.g. WordMorphism).
                         
    -  ``algorithm`` - string (default: ``'recursive'``) specifying which 
       algorithm to be used when computing the iterated palindromic closure. 
       It must be one of the two following values:

       - ``'definition'`` - computed using the definition
       - ``'recursive'`` - computation based on an efficient formula 
         that recursively computes the iterated right palindromic closure 
         without having to recompute the longest `f`-palindromic suffix 
         at each iteration [2].
                         
    OUTPUT:

        word -- the iterated (`f`-)palindromic closure of self
        
    EXAMPLES::

        sage: w = Word('abc')
        sage: w.iterated_right_palindromic_closure()
        word: abacaba

    ::

        sage: w = Word('aaa')
        sage: w.iterated_right_palindromic_closure()
        word: aaa

    ::

        sage: w = Word('abbab')
        sage: w.iterated_right_palindromic_closure()
        word: ababaabababaababa

    A right `f`-palindromic closure::

        sage: f = WordMorphism('a->b,b->a')
        sage: w = Word('abbab')
        sage: w.iterated_right_palindromic_closure(f=f)
        word: abbaabbaababbaabbaabbaababbaabbaab

    An infinite word::

        sage: t = words.ThueMorseWord('ab')
        sage: t.iterated_right_palindromic_closure()
        word: ababaabababaababaabababaababaabababaabab...

    There are two implementations computing the iterated right
    `f`-palindromic closure, the latter being much more efficient::

        sage: w = Word('abaab')
        sage: u = w.iterated_right_palindromic_closure(algorithm='definition')
        sage: v = w.iterated_right_palindromic_closure(algorithm='recursive')
        sage: u
        word: abaabaababaabaaba
        sage: u == v
        True
        sage: w = words.RandomWord(8)
        sage: u = w.iterated_right_palindromic_closure(algorithm='definition')
        sage: v = w.iterated_right_palindromic_closure(algorithm='recursive')
        sage: u == v
        True

    TESTS:

    The empty word::

        sage: w = Word()
        sage: w.iterated_right_palindromic_closure()
        word:

    If the word is finite, so is the result::

        sage: w = Word([0,1]*7)
        sage: c = w.iterated_right_palindromic_closure()
        sage: type(c)
        <class 'sage.combinat.words.word.FiniteWord_iter_with_caching'>

    REFERENCES:

    -   [1] A. de Luca, A. De Luca, Pseudopalindrome closure operators
        in free monoids, Theoret. Comput. Sci. 362 (2006) 282--300.      
    -   [2] J. Justin, Episturmian morphisms and a Galois theorem on
        continued fractions, RAIRO Theoret. Informatics Appl. 39 (2005)
        207-215.
    """
    from sage.combinat.words.word import FiniteWord_class, InfiniteWord_class
    if isinstance(self, FiniteWord_class):
        length = "finite"
    elif isinstance(self, InfiniteWord_class):
        length = None
    else:
        length = "unknown"
    if algorithm == 'definition':
        it = self._iterated_right_palindromic_closure_iterator(f=f)
    elif algorithm == 'recursive':
        it = self._iterated_right_palindromic_closure_recursive_iterator(f=f)
    else:
        raise ValueError, "algorithm (=%s) must be either 'definition' or\
                           'recursive'"
    return self._parent(it, length=length)
