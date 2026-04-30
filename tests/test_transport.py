import numpy as np

from gcascade_v5 import transport


def _diag():
    return transport.zero_diagnostics()


def test_gamma_pp_creates_electrons_and_conserves_nonnegativity():
    n = len(transport.energies)
    gamma = np.zeros(n)
    gamma[120] = 1.0
    electron = np.zeros(n)
    diagnostics = _diag()

    out_gamma, out_electron = transport.transport_step(
        gamma,
        electron,
        step_size=1.0,
        z_index=0,
        pp_log_extinction_row=np.full(n, np.log(0.5)),
        ics_log_extinction_row=np.zeros(n),
        d_edt_ics_row=np.zeros(n),
        pp_kernel=np.eye(n),
        pp_below_grid_row=np.zeros(n),
        ics_gamma_kernel=np.zeros((n, n)),
        ics_electron_kernel=np.eye(n),
        ics_gamma_energy_row=np.zeros(n),
        ics_below_grid_row=np.zeros(n),
        cel_data=None,
        diagnostics=diagnostics,
    )

    assert out_gamma[120] < gamma[120]
    assert out_electron[120] > 0.0
    assert np.all(out_gamma >= 0.0)
    assert np.all(out_electron >= 0.0)
    assert diagnostics["pp_absorbed_energy"] > 0.0


def test_electron_ics_creates_gamma_and_lower_energy_electrons():
    n = len(transport.energies)
    gamma = np.zeros(n)
    electron = np.zeros(n)
    electron[120] = 1.0
    diagnostics = _diag()
    shifted = np.zeros((n, n))
    shifted[119, 120] = 1.0

    out_gamma, out_electron = transport.transport_step(
        gamma,
        electron,
        step_size=1.0,
        z_index=0,
        pp_log_extinction_row=np.zeros(n),
        ics_log_extinction_row=np.full(n, np.log(0.5)),
        d_edt_ics_row=np.ones(n),
        pp_kernel=np.zeros((n, n)),
        pp_below_grid_row=np.zeros(n),
        ics_gamma_kernel=np.eye(n),
        ics_electron_kernel=shifted,
        ics_gamma_energy_row=np.full(n, transport.energies),
        ics_below_grid_row=np.zeros(n),
        cel_data=None,
        diagnostics=diagnostics,
    )

    assert out_gamma[120] > 0.0
    assert out_electron[119] > 0.0
    assert np.all(out_gamma >= 0.0)
    assert np.all(out_electron >= 0.0)
    assert diagnostics["ics_scattered_electron_energy"] > 0.0


def test_mixed_injection_is_additive_for_linear_step():
    n = len(transport.energies)
    gamma = np.zeros(n)
    electron = np.zeros(n)
    gamma[90] = 1.0
    electron[120] = 1.0

    kwargs = dict(
        step_size=1.0,
        z_index=0,
        pp_log_extinction_row=np.full(n, np.log(0.8)),
        ics_log_extinction_row=np.full(n, np.log(0.7)),
        d_edt_ics_row=np.ones(n),
        pp_kernel=np.eye(n),
        pp_below_grid_row=np.zeros(n),
        ics_gamma_kernel=np.eye(n),
        ics_electron_kernel=np.eye(n),
        ics_gamma_energy_row=np.full(n, transport.energies),
        ics_below_grid_row=np.zeros(n),
        cel_data=None,
    )
    g1, e1 = transport.transport_step(gamma, np.zeros(n), diagnostics=_diag(), **kwargs)
    g2, e2 = transport.transport_step(np.zeros(n), electron, diagnostics=_diag(), **kwargs)
    gb, eb = transport.transport_step(gamma, electron, diagnostics=_diag(), **kwargs)

    assert np.allclose(gb, g1 + g2)
    assert np.allclose(eb, e1 + e2)


def test_synchrotron_loss_is_recorded_without_negative_spectra():
    n = len(transport.energies)
    electron = np.zeros(n)
    electron[120] = 1.0
    diagnostics = _diag()

    out_gamma, out_electron = transport.transport_step(
        np.zeros(n),
        electron,
        step_size=1.0,
        z_index=0,
        pp_log_extinction_row=np.zeros(n),
        ics_log_extinction_row=np.full(n, np.log(0.5)),
        d_edt_ics_row=np.zeros(n),
        pp_kernel=np.zeros((n, n)),
        pp_below_grid_row=np.zeros(n),
        ics_gamma_kernel=np.zeros((n, n)),
        ics_electron_kernel=np.eye(n),
        ics_gamma_energy_row=np.zeros(n),
        ics_below_grid_row=np.zeros(n),
        cel_data=None,
        diagnostics=diagnostics,
        b_field_gauss=1.0e-7,
    )

    assert diagnostics["synchrotron_energy_lost"] > 0.0
    assert np.all(out_gamma >= 0.0)
    assert np.all(out_electron >= 0.0)


def test_synchrotron_step_conserves_energy_with_below_grid_term():
    n = len(transport.energies)
    electron = np.zeros(n)
    electron[220] = 1.0
    diagnostics = _diag()

    _, out_electron = transport.transport_step(
        np.zeros(n),
        electron,
        step_size=1.0,
        z_index=0,
        pp_log_extinction_row=np.zeros(n),
        ics_log_extinction_row=np.zeros(n),
        d_edt_ics_row=np.zeros(n),
        pp_kernel=np.zeros((n, n)),
        pp_below_grid_row=np.zeros(n),
        ics_gamma_kernel=np.zeros((n, n)),
        ics_electron_kernel=np.eye(n),
        ics_gamma_energy_row=np.zeros(n),
        ics_below_grid_row=np.zeros(n),
        cel_data=None,
        diagnostics=diagnostics,
        b_field_gauss=1.0e-2,
    )

    initial = transport.spectrum_energy(electron)
    final_known = (
        transport.spectrum_energy(out_electron)
        + diagnostics["synchrotron_energy_lost"]
        + diagnostics["below_grid_energy_lost"]
    )
    assert diagnostics["synchrotron_energy_lost"] > 0.0
    assert np.all(out_electron >= 0.0)
    assert abs(initial - final_known) <= 1.0e-10 * initial


def test_cel_transition_is_nonnegative_and_records_diagnostic():
    n = len(transport.energies)
    electron = np.zeros(n)
    electron[40] = 1.0
    mask = np.zeros(n, dtype=bool)
    mask[40] = True
    targets = transport.energies.copy()
    targets[40] = transport.energies[39]
    diagnostics = _diag()

    _, out_electron = transport.transport_step(
        np.zeros(n),
        electron,
        step_size=1.0,
        z_index=0,
        pp_log_extinction_row=np.zeros(n),
        ics_log_extinction_row=np.full(n, np.log(0.5)),
        d_edt_ics_row=np.ones(n),
        pp_kernel=np.zeros((n, n)),
        pp_below_grid_row=np.zeros(n),
        ics_gamma_kernel=np.zeros((n, n)),
        ics_electron_kernel=np.zeros((n, n)),
        ics_gamma_energy_row=np.zeros(n),
        ics_below_grid_row=np.zeros(n),
        cel_data=(mask, targets),
        diagnostics=diagnostics,
    )

    assert np.all(out_electron >= 0.0)
    assert out_electron[39] > 0.0
    assert diagnostics["cel_scattered_electron_energy"] > 0.0
