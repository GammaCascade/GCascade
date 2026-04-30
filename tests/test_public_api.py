import pytest

import gcascade_v5 as g


def test_public_functions_exist():
    for name in [
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
        "CascadeResult",
        "get_bundle_info",
        "reset_factory_settings",
    ]:
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
    assert g.change_ebl_model is g.changeEBLModel
    assert g.change_magnetic_field is g.changeMagneticField


def test_path_helpers_and_setters(tmp_path):
    lib = (tmp_path / "LibrariesV5").resolve()
    lib.mkdir(parents=True, exist_ok=True)

    g.set_library_path(lib)
    assert g.get_library_path() == lib

    g.reset_state()


def test_set_library_path_rejects_missing(tmp_path):
    with pytest.raises(FileNotFoundError):
        g.set_library_path(tmp_path / "missing_libraries_v5")


def test_public_version_string():
    assert g.__version__ == "5.1"


def test_eblindex_tracks_change(quick_bundle_path):
    g.set_library_path(quick_bundle_path)
    original = g.EBLindex
    try:
        g.changeEBLModel(6)
        assert g.EBLindex == 6
    finally:
        if g.EBLindex != original:
            g.changeEBLModel(original)

