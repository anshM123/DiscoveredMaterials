"""Physics invariants of the framework (run: python -m pytest -q tests/).
Cheap: uses cached phonons only."""
import os, sys
import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from phonons import load, DW, W, channel_spectra          # noqa: E402
from model import Side, Interface, G_DMM_debye, G_AMM_debye, coincidence_cell   # noqa: E402
from predict import load_params, G_ep                      # noqa: E402

CRYS = [("Si", "111"), ("GaN", "0001"), ("Al", "111"), ("AlN", "0001"), ("SiC4H", "0001"), ("Al2O3", "0001")]
P = load_params()


@pytest.mark.parametrize("name,normal", CRYS)
def test_dos_normalisation(name, normal):
    d = load(name)
    assert abs(np.sum(d["D"]) * DW * d["V"] / (3 * d["nat"]) - 1) < 1e-3


@pytest.mark.parametrize("name,normal", [("Si", "111"), ("GaN", "0001"), ("Si", "100")])
def test_identical_crystals_transmit_fully(name, normal):
    s = Side(name, normal)
    I = Interface(s, Side(name, normal))
    for T in (100, 300, 700):
        assert abs(I.G_elastic(T) / s.G_rad(T) - 1) < 2e-3          # no interface -> radiation limit
        assert abs(I.G_pp(T, P["alpha"], P["b"]) / s.G_rad(T) - 1) < 2e-3


@pytest.mark.parametrize("a,b", [(("AlN", "0001"), ("GaN", "0001")), (("Al", "111"), ("Si", "111")),
                                 (("GaAs", "100"), ("AlAs", "100")), (("Au", None), ("PMMA", None))])
def test_reciprocity_and_bounds(a, b):
    s1, s2 = Side(*a), Side(*b)
    I12, I21 = Interface(s1, s2), Interface(s2, s1)
    assert np.all(I12.tau_el >= -1e-12) and np.all(I12.tau_el <= 1 + 1e-9)
    for T in (100, 300, 800):
        g12, g21 = I12.G_elastic(T), I21.G_elastic(T)
        assert abs(g12 / g21 - 1) < 1e-6                                # detailed balance
        assert g12 <= min(s1.G_rad(T), s2.G_rad(T)) * (1 + 1e-6)       # Landauer bound
        assert I12.G_pp(T, P["alpha"], P["b"]) >= g12                   # inelastic channel only adds


@pytest.mark.parametrize("name,normal", [("GaN", "0001"), ("Si", "111")])
def test_channel_spectra_sum_to_flux(name, normal):
    Mk = channel_spectra(name, normal, nk=16, nz=32)
    assert abs(Mk.mean(axis=0).sum() / load(name)[f"M_{normal}"].sum() - 1) < 0.07


def test_coincidence_known_cases():
    assert abs(coincidence_cell(Side("Al", "111"), Side("Si", "111"))[0] - 11.5) < 0.2      # 4:3
    assert abs(coincidence_cell(Side("AlN", "0001"), Side("GaN", "0001"))[0] - 3.2) < 0.1   # 1:1
    assert coincidence_cell(Side("Au", "111"), Side("PMMA")) is None                         # amorphous


def test_baselines_positive_and_ordered():
    s1, s2 = Side("Al", "111"), Side("Si", "111")
    for T in (100, 300):
        assert G_DMM_debye(s1, s2, T) > 0 and G_AMM_debye(s1, s2, T) > 0


def test_eph_monotonic_and_missing_metal():
    assert G_ep("Al", 300) > G_ep("Al", 700) > 0
    with pytest.warns(UserWarning):
        assert G_ep("Unobtainium", 300) == np.inf


def test_agf_blocks_hermitian():
    pytest.importorskip("phonopy")
    from phonopy import load as phload
    import agf
    path = os.path.join(agf.CACHE, "GaAs100_AlAs100_ph.yaml")
    if not os.path.exists(path):
        pytest.skip("AGF force constants not cached")
    uL, uR = agf.oriented_unit("GaAs", "100"), agf.oriented_unit("AlAs", "100")
    pl = [max(1, int(np.ceil(agf.PL_MIN / u.cell[2, 2]))) for u in (uL, uR)]
    sl, tag, _ = agf.superlattice(uL, 4 * pl[0], uR, 4 * pl[1])
    A = agf.AGF(phload(path, log_level=0), tag, 4, tuple(pl))
    for k in agf.kmesh(A.rec2d, 4)[:3]:
        H00L, H01L, H00R, H01R, Hd, VL, VR, sizes = A.lead_and_device(k)
        for H in (H00L, H00R, Hd):
            assert np.abs(H - H.conj().T).max() < 1e-12
