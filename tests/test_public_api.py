import pytest

import gcascade_v5 as g
import gcascade_v5.builders as builders


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
        "convert_legacy_library",
        "get_bundle_info",
        "list_generated_variants",
        "get_active_cycle_path",
        "set_active_cycle_path",
        "reset_factory_settings",
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
    assert g.change_ebl_model is g.changeEBLModel
    assert g.change_magnetic_field is g.changeMagneticField


def test_legacy_namespace_is_exposed():
    assert hasattr(g, "legacy")
    assert hasattr(g.legacy, "CascadePoint")


def test_path_helpers_and_setters(tmp_path):
    lib = (tmp_path / "LibrariesV5").resolve()
    generated = (tmp_path / "generated_libraries").resolve()
    lib.mkdir(parents=True, exist_ok=True)

    g.set_library_path(lib)
    g.set_generated_library_path(generated)

    assert g.get_library_path() == lib
    assert g.get_generated_library_path() == generated

    g.reset_state()


def test_set_library_path_rejects_missing(tmp_path):
    with pytest.raises(FileNotFoundError):
        g.set_library_path(tmp_path / "missing_libraries_v5")


def test_numba_status_helpers_are_consistent():
    assert isinstance(g.is_numba_available(), bool)
    assert isinstance(g.get_numba_enabled(), bool)
    if not g.is_numba_available():
        assert g.get_numba_enabled() is False


def test_public_version_string():
    assert g.__version__ == "5.0"


def test_generated_variant_filename_format():
    assert builders._variant_id("Dom", 1.0e-7, 0.0) == "ebl_Dom_B_1e-7_gamma_0e0"
    assert builders._variant_id("SL", 1.2345e-7, 1.2345) == "ebl_SL_B_1.23e-7_gamma_1.23e0"


def test_eblindex_tracks_change():
    original = g.EBLindex
    try:
        g.changeEBLModel(6)
        assert g.EBLindex == 6
    finally:
        if g.EBLindex != original:
            g.changeEBLModel(original)
