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

########################################################################
# Functions related to the computation of iterated palindromic closure #
########################################################################

def occurrences_iterator(self, other, reverse=False):
    r"""
    Returns an iterator enumerating all occurrences of
    other in self.

    INPUT:

    - ``other`` - the pattern that needs to be matched.
    - ``reverse`` - boolean (default: False). When set to True,
      the search is performed from the right to the left.

    OUTPUT:

    iterator

    EXAMPLES::

        sage: w = words.FibonacciWord()[:20]
        sage: w
        word: 01001010010010100101
        sage: w[:3]
        word: 010
        sage: it = occurrences_iterator(w, w[:3])
        sage: list(it)
        [0, 3, 5, 8, 11, 13, 16]
        sage: it = occurrences_iterator(w, w[:3], reverse=True)
        sage: list(it)
        [16, 13, 11, 8, 5, 3, 0]

    ::

        sage: w[1:2]
        word: 1
        sage: it = occurrences_iterator(w, w[1:2])
        sage: list(it)
        [1, 4, 6, 9, 12, 14, 17, 19]
        sage: it = occurrences_iterator(w, w[1:2], reverse=True)
        sage: list(it)
        [19, 17, 14, 12, 9, 6, 4, 1]
    """
    n = len(self)
    if reverse:
        i = n - 1
    else:
        i = 0
    while (not reverse and i < n) or (reverse and i >=0):
        p = len(other)
        if p <= n - i:
            if self[i:i+p] == other:
                yield i
        if reverse:
            i -= 1
        else:
            i += 1
    raise StopIteration

def _rightmost_occurrence_preceded_by_closure_type(self, involutions, letter, involution):
    r"""
    Returns the rightmost occurrence of the given letter
    in self preceded by the given involution in involutions.
    Returns None on failure.

    INPUT:

    - ``involutions`` - iterable enumerating the successive
      involutions that needs to be applied for computing the
      iterated right palindromic closure.
    - ``letter`` - the letter to be searched for.
    - ``involution`` - the involution that we want the letter
      to be preceded by.
    """
    occ_iter = occurrences_iterator(self, letter, reverse=True)
    for o in occ_iter:
        if o > 0 and involutions[o - 1] == involution:
            return o
    return None

def _iterated_right_palindromic_closure_naive_iterator(self, involutions):
    r"""
    Returns an iterator over the right palindromic closure
    of self with respect to the given involutions sequence
    using the naive algorithm.
    """
    parent = self.parent()
    w = self[:0]
    for (letter, involution) in zip(self, involutions):
        length_before = w.length()
        w = (w * parent([letter])).palindromic_closure(f=involution)
        length_after = w.length()
        d = length_after - length_before
        for a in w[-d:]:
            yield a
    raise StopIteration

def _iterated_right_palindromic_closure_recursive_iterator(self, involutions):
    r"""
    Returns an iterator over the right palindromic closure
    of self with respect to the given involutions sequence
    using an efficient formula.
    The sequence self must be binary and it is assumed to be
    normalized as well.
    """
    parent = self.parent()
    if len(parent.alphabet()) != 2:
        raise ValueError, 'self (=%s) must be over a binary alphabet'
    a,b = parent.alphabet()
    id = WordMorphism({a : a, b : b}, codomain=parent)
    ex = WordMorphism({a : b, b : a}, codomain=parent)
    w = self[:0]
    lengths = [0]
    for (i, letter, involution) in zip(range(len(self)), self, involutions):
        if i == 0 or i == 1:
            w += parent([letter])
        else:
            o = _rightmost_occurrence_preceded_by_closure_type(self[:i], involutions[:i], involution(involutions[i - 1](letter)), involution)
            if involution == id and involutions[i - 1] == id:
                if letter not in self[:i]:
                    w += parent([letter]) + w
                elif o is not None:
                    w += w[lengths[o]:]
                else:
                    w += w
            elif involution == id and involutions[i - 1] == ex:
                print i, o
                if o is not None:
                    w += ex(w[lengths[o]:])
                else:
                    w += ex(w)
            elif involution == ex and involutions[i - 1] == id:
                if letter not in self[:i]:
                    w += ex(w)
                elif o is not None:
                    w += ex(w[lengths[o]:])
                else:
                    w += parent([letter]) + ex(parent([letter])) + ex(w)
            elif involution == ex and involutions[i - 1] == ex:
                if o is not None:
                    w += w[lengths[o]:]
                else:
                    w += parent([letter]) + ex(parent([letter])) + w
            print w
        lengths.append(len(w))
        for a in w[lengths[i]:]:
            yield a
    raise StopIteration

def _iterated_right_palindromic_closure_recursive_iterator_merged(self, involutions):
    r"""
    Returns an iterator over the right palindromic closure
    of self with respect to the given involutions sequence
    using an efficient formula.
    The sequence self must be binary and it is assumed to be
    normalized as well.
    """
    parent = self.parent()
    if len(parent.alphabet()) != 2:
        raise ValueError, 'self (=%s) must be over a binary alphabet'
    a,b = parent.alphabet()
    id = WordMorphism({a : a, b : b}, codomain=parent)
    ex = WordMorphism({a : b, b : a}, codomain=parent)
    w = self[:0]
    lengths = [0]
    for (i, letter, involution) in zip(range(len(self)), self, involutions):
        if letter not in self[:i] or involution not in involutions[:i]:
            if involution == id:
                w += parent([letter]) + w
            elif involution == ex:
                if letter == w[-1]:
                    w += parent([letter]) + ex(parent([letter])) + ex(w)
                else:
                    w += ex(w)
        else:
            o = _rightmost_occurrence_preceded_by_closure_type(self[:i], involutions[:i], involution(involutions[i - 1](letter)), involution)
            if o is None:
                if involution == id:
                    w += involution(involutions[i - 1](w))
                elif involution == ex:
                    if letter == w[-1]:
                        w += parent([letter]) + ex(parent([letter])) + involution(involutions[i - 1](w))
                    else:
                        w += involution(involutions[i - 1](w))

            else:
                w += involution(involutions[i - 1](w[lengths[o]:]))
        lengths.append(len(w))
        for a in w[lengths[i]:]:
            yield a
    raise StopIteration

def iterated_right_palindromic_closure(self, involutions=None, algorithm='naive'):
    r"""
    Returns the iterated palindromic closure of self with
    respect to the given involutions sequence.

    INPUT:

    - ``involutions`` - iterable (default: None) on involutions

    OUTPUT:

    word

    EXAMPLE::

        sage: A = Words([0,1])
        sage: R = WordMorphism({0:0,1:1}, codomain=A)
        sage: E = WordMorphism({0:1,1:0}, codomain=A)
        sage: T = [R,E,E,E,E,E,R,R,R,E] + [R,E]
        sage: w = A([0,0,0,1,1,0,1,0,1,1] + [0,0])
        sage: u = iterated_right_palindromic_closure(w, T)
        sage: v = u[:2000]
        sage: [round(float(v.number_of_factors(i) / (4 * i)),3) for i in range(1,150)]
        [0.5, 0.5, 0.5, 0.5, 0.5, 0.583, 0.643, 0.688, 0.722, 0.8, 0.864, 0.917, 0.962, 1.0, 1.033, 1.063, 1.088, 1.111, 1.132, 1.15, 1.167, 1.182, 1.196, 1.208, 1.22, 1.212, 1.204, 1.196, 1.19, 1.183, 1.177, 1.172, 1.167, 1.162, 1.157, 1.153, 1.149, 1.145, 1.141, 1.138, 1.134, 1.131, 1.128, 1.125, 1.122, 1.12, 1.117, 1.115, 1.112, 1.11, 1.108, 1.106, 1.104, 1.102, 1.1, 1.098, 1.096, 1.095, 1.093, 1.092, 1.09, 1.089, 1.087, 1.086, 1.085, 1.083, 1.082, 1.081, 1.08, 1.079, 1.077, 1.076, 1.075, 1.074, 1.073, 1.072, 1.071, 1.071, 1.07, 1.069, 1.068, 1.067, 1.066, 1.065, 1.065, 1.064, 1.063, 1.063, 1.062, 1.061, 1.06, 1.06, 1.059, 1.059, 1.058, 1.057, 1.057, 1.056, 1.056, 1.055, 1.054, 1.054, 1.053, 1.053, 1.052, 1.052, 1.051, 1.051, 1.05, 1.05, 1.05, 1.049, 1.049, 1.048, 1.048, 1.047, 1.047, 1.047, 1.046, 1.046, 1.045, 1.045, 1.045, 1.044, 1.044, 1.044, 1.043, 1.043, 1.043, 1.042, 1.042, 1.042, 1.041, 1.041, 1.041, 1.04, 1.04, 1.04, 1.04, 1.039, 1.039, 1.039, 1.038, 1.038, 1.038, 1.038, 1.037, 1.037, 1.037]
    """
    if involutions is None:
        alphabet = self.parent().alphabet()
        id = WordMorphism(lambda a: a)
        from itertools import cycle
        involutions = cycle(id)
    if algorithm == 'naive':
        iter = _iterated_right_palindromic_closure_naive_iterator(self, involutions)
    elif algorithm == 'recursive':
        iter = _iterated_right_palindromic_closure_recursive_iterator(self, involutions)
    elif algorithm == 'merged':
        iter = _iterated_right_palindromic_closure_recursive_iterator_merged(self, involutions)
    else:
        raise ValueError, "algorithm parameter (=%s) must be 'naive' or 'recursive'"
    return Word(iter)


##########################################
# Functions related to the normalization #
##########################################

def iterated_right_palindromic_closure_lengths_iterator(self, involutions):
    r"""
    Returns an iterator over the length of the successive
    palindromes obtained when computing the right
    palindromic closure of self with respect to the given
    involutions sequence using the naive algorithm.
    """
    parent = self.parent()
    w = self[:0]
    for (letter, involution) in zip(self, involutions):
        length_before = w.length()
        yield length_before
        w = (w * parent([letter])).palindromic_closure(f=involution)
        length_after = w.length()
        d = length_after - length_before
    yield len(w)

def palindrome_prefixes_lenghts(self):
    r"""
    Returns the lengths of the f-palindrome prefixes of
    self for all possible f
    """
    alphabet = set(self)
    if len(alphabet) > 2:
        raise ValueError, 'self (=%s) must be over a binary alphabet'
    elif len(alphabet) == 1:
        return range(len(self) + 1)
    a,b = alphabet
    ex = WordMorphism({a : b, b : a})
    lengths = []
    for prefix in self.prefixes_iterator():
        if prefix.is_palindrome() or prefix.is_palindrome(f=ex):
            lengths.append(len(prefix))
    return lengths

def is_normalized(self, involutions):
    r"""
    Returns True if self is normalized, i.e. that when
    computing its iterated right palindromic closure,
    its f-palindrome prefixes are exactly those given
    by the closure.
    """
    return len(list(iterated_right_palindromic_closure_lengths_iterator(self, involutions))) == len(palindrome_prefixes_lenghts(iterated_right_palindromic_closure(self, involutions)))

def is_normalized_fast(self, involutions):
    r"""
    Returns True if self is normalized by looking at
    the form of the word
    """
    pass

def not_normalized_words_iterator(max_length=None):
    r"""
    Returns an iterator over the not normalized words
    """
    if max_length is None:
        max_length = +Infinity
    alphabet = [0, 1]
    parent = Words(alphabet)
    a,b = alphabet
    id = WordMorphism({a : a, b : b}, codomain=parent)
    ex = WordMorphism({a : b, b : a}, codomain=parent)
    from heapq import heappop, heappush
    queue = []
    if max_length >= 1:
        for morphism in [id, ex]:
            queue.append((1, parent([a]), [morphism]))
    while queue:
        (i,w,m) = heappop(queue)
        if not is_normalized(w, m):
            yield (w, m)
        elif len(w) < max_length:
            for letter in [a, b]:
                for morphism in [id, ex]:
                    heappush(queue, (i + 1, w + parent([letter]), m + [morphism]))
    raise StopIteration

def normalized_words_iterator(max_length=None):
    r"""
    Returns an iterator over the normalized words
    """
    if max_length is None:
        max_length = +Infinity
    alphabet = [0, 1]
    parent = Words(alphabet)
    a,b = alphabet
    id = WordMorphism({a : a, b : b}, codomain=parent)
    ex = WordMorphism({a : b, b : a}, codomain=parent)
    from heapq import heappop, heappush
    queue = []
    if max_length >= 1:
        for morphism in [id, ex]:
            queue.append((1, parent([a]), [morphism]))
    while queue:
        (i,w,m) = heappop(queue)
        if is_normalized(w, m):
            yield (w, m)
            if len(w) < max_length:
                for letter in [a, b]:
                    for morphism in [id, ex]:
                        heappush(queue, (i + 1, w + parent([letter]), m + [morphism]))
    raise StopIteration

def display_nice_list_of_not_normalized_words(max_length=None):
    r"""
    Displays a list of not normalized words
    """
    alphabet = [0, 1]
    parent = Words(alphabet)
    a,b = alphabet
    id = WordMorphism({a : a, b : b}, codomain=parent)
    ex = WordMorphism({a : b, b : a}, codomain=parent)
    it = not_normalized_words_iterator(max_length)
    for (w,m) in it:
        s = ''
        for letter, morphism in zip(w, m):
            if morphism(a)[0] != a:
                s += '^'
            s += str(letter)
        print s    

##############################
# Functions to detect period #
##############################

def minimal_period(self, f=None):
    r"""
    Returns the smallest `f`-period of self.

    EXAMPLES::

        sage: ex = WordMorphism('a->b,b->a')
        sage: w  = Word('ababababa')
        sage: minimal_period(w)
        2
        sage: minimal_period(w, ex)
        1
    """
    if f is None:
        f = lambda x:x
    p = len(self)
    for i in range(p, 0, -1):
        even = False
        j = i
        prefix = self[:i]
        has_period_i = True
        while has_period_i and j < len(self):
            if j + i < len(self):
                w = self[j:j+i]
            else:
                w = self[j:]
            if not even:
                w = f(w)
            even = not even
            has_period_i = w.is_prefix(prefix)
            j += i
        if has_period_i:
            p = i
    return p

###############################
# Functions to test algorithm #
###############################

def test_algorithm(max_length):
    it = normalized_words_iterator(max_length)
    identify = lambda i:'id' if i.is_identity() else 'ex'
    for (w, involutions) in it:
        merged = iterated_right_palindromic_closure(w, involutions, algorithm='merged')
        naive = iterated_right_palindromic_closure(w, involutions, algorithm='naive')
        if merged != naive:
            print 'problem with ' + str((w,map(identify, involutions)))
            print 'naive  = %s'%naive
            print 'merged = %s'%merged

def test_bisequence(u, v):
    alphabet = [0, 1]
    parent = Words(alphabet)
    a,b = alphabet
    id = WordMorphism({a : a, b : b}, codomain=parent)
    ex = WordMorphism({a : b, b : a}, codomain=parent)
    w = parent(map(lambda L:int(L), u))
    involutions = map(lambda m: id if m == 'R' else ex, v)
    return iterated_right_palindromic_closure(w, involutions, algorithm='merged')

#######################
# Bezout coefficients #
#######################

def bezout(a, b):
    r = a; rp = b
    u = 1; up = 0
    v = 0; vp = 1
    while rp != 0:
        q = r // rp
        rs = r; us = u; vs = v
        r = rp; u = up; v = vp
        rp = rs - q * rp
        up = us - q * up
        vp = vs - q * vp
    return (u, v)

def enumerate_impossible_cases(n):
    for i in range(n):
        for j in range(i + 1, n):
            (x,y) = bezout(i, j) 
            if (x + y) % 2 == 0:
                print i, j, gcd(i, j)

####################
# Class bisequence #
####################

class Bisequence:

    #################
    # Basic methods #
    #################

    def __init__(self, directive_word, closure_type_word):
        r"""
        Creates a directive bisequence of a generalized
        pseudostandard word.

        INPUT:

        - ``directive_word`` - binary string on `\{0,1\}`
        - ``closure_type_word`` - binary string on `\{R,E\}`

        EXAMPLES::

            sage: s = Bisequence('011', 'RER')
            sage: s
            bisequence: (011, RER)
        """
        if type(directive_word) != str or type(closure_type_word) != str:
            raise ValueError, 'directive word (=%s) and closure_type_word (=%s) must be strings.'
        elif set(directive_word) != set(['0', '1']):
            raise ValueError, 'directive word (=%s) must be on alphabet {0,1}.'
        elif set(closure_type_word) != set(['R', 'E']):
            raise ValueError, 'closure type word (=%s) must be on alphabet {R,E}.'
        self._directive_word = directive_word
        self._closure_type_word = map(self._involutions_map, closure_type_word)

    def _repr_(self):
        r"""
        Returns a string representation of self.

        TEST::

            sage: Bisequence('011', 'RER')
            bisequence: (011, RER)
        """
        return 'bisequence: (%s, %s)'%(self._directive_word, self._closure_type_word)

    __repr__ = _repr_

    def directive_word(self):
        r"""
        Returns the directive word of self.
        """
        return self._directive_word

    def closure_type_word(self):
        r"""
        Returns the closure type word of self.
        """
        return self._closure_type_word

    #################
    # Normalization #
    #################

    def _iterated_right_palindromic_closure_lengths_iterator(self):
        r"""
        Returns an iterator over the length of the successive
        palindromes obtained when computing the right
        palindromic closure of self.
        """
        dw = self.directive_word()
        ctw = self.closure_type_word()
        parent = dw.parent()
        w = directive_word[:0]
        for (letter, involution) in zip(dw, ctw):
            length_before = w.length()
            yield length_before
            w = (w * parent([letter])).palindromic_closure(f=involution)
            length_after = w.length()
            d = length_after - length_before
        yield len(w)

    def is_normalized(self, algorithm='naive'):
        r"""
        Returns True is self is a normalized bisequence,
        i.e. if it does not miss any palindrome or
        pseudopalindrome when computing the iterated
        right palindromic closure.
        """
        pass
