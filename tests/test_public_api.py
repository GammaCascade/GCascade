import pytest

import gcascade_v5 as g


def test_public_functions_exist():
    names = [
        "RedshiftPoint",
        "AttenuatePoint",
        "CascadePoint",
        "RedshiftDiffuse",
        "AttenuateDiffuse",
        "CascadeDiffuse",
        "RedshiftEvolving",
        "AttenuateEvolving",
        "CascadeEvolving",
        "changeEBLModel",
        "changeMagneticField",
    ]
    for name in names:
        assert hasattr(g, name)


def test_snake_case_aliases_exist():
    assert g.redshift_point is g.RedshiftPoint
    assert g.attenuate_point is g.AttenuatePoint
    assert g.cascade_point is g.CascadePoint
    assert g.redshift_diffuse is g.RedshiftDiffuse
    assert g.attenuate_diffuse is g.AttenuateDiffuse
    assert g.cascade_diffuse is g.CascadeDiffuse
    assert g.redshift_evolving is g.RedshiftEvolving
    assert g.attenuate_evolving is g.AttenuateEvolving
    assert g.cascade_evolving is g.CascadeEvolving


def test_path_helpers_and_setters(tmp_path):
    lib = (tmp_path / "LibrariesV4").resolve()
    generated = (tmp_path / "generated_libraries").resolve()
    lib.mkdir(parents=True, exist_ok=True)

    g.set_library_path(lib)
    g.set_generated_library_path(generated)

    assert g.get_library_path() == lib
    assert g.get_generated_library_path() == generated

    g.reset_state()


def test_set_library_path_rejects_missing(tmp_path):
    with pytest.raises(FileNotFoundError):
        g.set_library_path(tmp_path / "missing_libraries_v4")
