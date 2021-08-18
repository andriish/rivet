import pytest
from old_plotfile_converter import parse_old_plotfile, type_conversion


@pytest.mark.parametrize(
    's,expected,expected_type',
    [
        ('singleWord', 'singleWord', str), 
        ('multiple words', 'multiple words', str), 
        ('string and numbers 123', 'string and numbers 123', str), 
        (r'Pseudorapidity $\sqrt{s}=7$ TeV, INEL $>0$', r'Pseudorapidity $\sqrt{s}=7$ TeV, INEL $>0$', str), 
        ('1', 1, int), 
        ('0', 0, int),
        ('87', 87, int),
        ('-10', -10, int),
        ('10.', 10., float),
        ('[1]', [1], list),
        ('[1, 2, 3]', [1, 2, 3], list),
    ]
)
def test_type_conversion(s, expected, expected_type):
    converted = type_conversion(s)
    assert converted == expected
    assert isinstance(converted, expected_type)


@pytest.mark.parametrize('s', 
    [0, 1, 87, -1, (1, 2, 3), [1, 2, 3], {'a':1}, False, True, None]
)
def test_type_conversion_error(s):
    with pytest.raises(TypeError):
        type_conversion(s)
