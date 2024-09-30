import pytest

from lsme.numeric.matchingresult import GenericMatchingResult
from lsme.numeric import GranadazetaMatchingResult as Zeta

def test_generic_matching_result_initialization():
    # Create an instance of GenericMatchingResult
    result = GenericMatchingResult(name="TestModel")

    # Check name is set correctly
    assert result.name == "TestModel"

    # Check default scale is 1.0
    assert result.scale == 1.0

    # Check e.g. that vev is present
    assert result.vev

def test_coeff_dict():
    import wilson
    import flavio
    from wilson import wcxf

    zeta = Zeta()
    coeff_dict_ = zeta.coeff_dict()

    my_wilson = wilson.Wilson(coeff_dict_, scale=1e3, eft='SMEFT', basis='Warsaw')
    # Validate the whole coeff dict
    my_wilson.wc.validate()

    # Replace just with one value for testing
    dict_ = {"lq3_1111": coeff_dict_["lq3_1111"]}
    my_wilson = wilson.Wilson(dict_, scale=1e3, eft='SMEFT', basis='Warsaw')

    # Run down to 160 GeV and fish out second-gen operator coeff which should
    # now be non-zero
    wc = my_wilson.match_run(scale=160, eft='SMEFT', basis='Warsaw')
    low_scale = wc["lq3_2222"]

    assert abs(low_scale.real) > 0
