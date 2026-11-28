"""Offline tests: no network, no seismic tools. `pytest -q tests`."""
import io
import os
import sys
import types

import numpy as np
import pandas as pd
import pytest
from obspy import Stream, Trace, UTCDateTime

import pqglobal
from pqglobal import catalog, network_picks, run, stations, velmodel, waveforms

COMCAT = """time,latitude,longitude,depth,mag,magType,id
2015-10-19T23:21:17.680Z,37.7918333,-121.9635,7.4,3.6,mw,nc72540765
2015-10-19T23:21:17.990Z,37.79,-121.96,7.0,1.1,md,nc72540766
2015-10-19T23:21:11.380Z,37.7905,-121.9640,7.3,2.32,md,nc72540760
"""


def test_package_finds_pocketquake_and_pipeline():
    import pocketquake
    assert os.path.isdir(os.path.join(pqglobal.EQCYCLE_DIR, "pipeline"))
    from pipeline.clusters._base import global_cluster, layered_velmodel  # noqa: F401
    assert pocketquake.__version__ >= "1.16.0"


def test_comcat_csv_and_same_second_rule(tmp_path):
    df = catalog.from_usgs_csv(io.StringIO(COMCAT))
    assert list(df.columns) == catalog.COLUMNS and len(df) == 3
    out = tmp_path / "event_catalog.csv"
    assert catalog.write(df, str(out)) == {"kept": 2, "dropped": 1}
    kept = pd.read_csv(out)
    assert set(kept.source_id) == {"nc72540765", "nc72540760"}          # the M1.1 in the M3.6's second goes
    assert list(kept.source_id) == ["nc72540760", "nc72540765"]          # time order
    assert pd.read_csv(str(out)[:-4] + ".dropped.csv").source_id.tolist() == ["nc72540766"]


def test_event_ids_match_the_pipeline_rule():
    assert waveforms.event_id("2015-10-19T23:21:17.680Z") == "20151019232117"


def test_pick_service_and_id_numbers():
    cat = pd.DataFrame({"source_id": ["nc1", "nc2", "ci3"]})
    assert network_picks.service_for(cat) == "NCEDC"
    assert network_picks._num("quakeml:nc.anss.org/Event/NC/72540765") == "72540765"
    assert network_picks._num("nc72540765") == "72540765"


def test_sensor_precedence():
    bands = ["HN", "EH", "HH", "CN", "EP"]
    assert sorted(bands, key=stations._rank) == ["HH", "EH", "EP", "HN", "CN"]
    assert stations.is_accelerometer("HN") and not stations.is_accelerometer("EP")


def test_con_cut_is_20_layers_within_30_ms():
    tops, vp = velmodel.cut(velmodel.CON_NODES, (5, 8, 5, 1))
    assert len(tops) == velmodel.MAX_LAYERS and tops[0] == 0.0 and vp[-1] == 7.95
    depths, dists = np.arange(1.0, 15.01, 1.0), np.arange(0.0, 61.0, 4.0)
    ref = velmodel.reference_times(velmodel.CON_NODES, depths, dists)
    err = velmodel.errors(tops, vp, depths, dists, ref)
    assert np.abs(err).max() < 30.0                                    # ms, P


def test_velocity_model_argument_builds_the_model():
    from pipeline.clusters._base import layered_velmodel
    m = eval(run._velocity_model_arg("con", (5, 8, 5, 1), 1.747), {"layered_velmodel": layered_velmodel})
    assert m.name == "con" and len(m.p_rows) == 20
    assert m.s_rows[0][0] == pytest.approx(m.p_rows[0][0] / 1.747, abs=1e-4)


@pytest.mark.parametrize("extra, expected", [
    ([], dict(s_on_t=False, orientation_csv=None, orientation_missing="error",
              rotation_reference=None)),
    (["--s-on-t"], dict(s_on_t=True, orientation_csv=None, orientation_missing="error",
                        rotation_reference=None)),
    (["--orientation-csv", "residual.csv", "--orientation-missing", "skip",
      "--rotation-reference", "37.79,-121.96"],
     dict(s_on_t=True, orientation_csv="residual.csv", orientation_missing="skip",
          rotation_reference=(37.79, -121.96))),
])
def test_cli_forwards_orientation_options(tmp_path, monkeypatch, extra, expected):
    """Cached FDSN stages must carry the requested orientation policy into PocketQuake."""
    catalog.write(catalog.from_usgs_csv(io.StringIO(COMCAT)), str(tmp_path / "event_catalog.csv"))
    (tmp_path / "network_picks.csv").touch()
    st_dir = tmp_path / "station_table"
    st_dir.mkdir()
    for name in run.STATION_FILES:
        (st_dir / name).touch()
    pd.DataFrame(dict(verdict=["used"], band=["HH"], location=[""])).to_csv(
        st_dir / "sensors.csv", index=False)
    called = []
    orchestrator = types.ModuleType("pocketquake.orchestrate")

    def capture(*args, **kwargs):
        called.append((args, kwargs))
        return {"status": "captured"}

    orchestrator.orchestrate = capture
    monkeypatch.setitem(sys.modules, "pocketquake.orchestrate", orchestrator)
    # Every FDSN stage is cached; any accidental network fetch is a test failure.
    def unexpected_fetch(*args, **kwargs):
        raise AssertionError("cached CLI test attempted a network fetch")

    monkeypatch.setattr(run.catmod, "from_usgs", unexpected_fetch)
    monkeypatch.setattr(run.network_picks, "fetch", unexpected_fetch)
    monkeypatch.setattr(run.stations, "select_picked", unexpected_fetch)
    run.main(["--name", "offline", "--region", "Offline", "--work", str(tmp_path),
              "--catalog-csv", "unused.csv", "--velmodel", "con", "--apply", *extra])
    assert len(called) == 1
    assert {key: called[0][1][key] for key in expected} == expected


@pytest.mark.parametrize("value", ["37.79", "37.79,-121.96,0", "nan,0", "91,0", "0,181"])
def test_rotation_reference_rejects_invalid_coordinates(value):
    import argparse
    with pytest.raises(argparse.ArgumentTypeError):
        run._rotation_reference(value)


def _tr(ch, t0, n, sr=100.0):
    return Trace(np.arange(n, dtype=float), header=dict(network="NC", station="X", channel=ch,
                                                        starttime=t0, sampling_rate=sr))


def test_common_stretch_keeps_the_origin_segment():
    o = UTCDateTime("2015-10-19T23:21:17")
    st = Stream([_tr("HHZ", o - 240, 60000), _tr("HHN", o - 240, 60000),
                 _tr("HHE", o - 240, 20000), _tr("HHE", o - 20, 38000)])   # E: gap, then a piece over o
    out = waveforms._common_stretch(st, o)
    assert out is not None and len(out) == 3
    assert len({t.stats.npts for t in out}) == 1
    assert all(t.stats.starttime <= o <= t.stats.endtime for t in out)
    assert out[0].stats.starttime == o - 20


def test_common_stretch_refuses_two_components():
    o = UTCDateTime(0)
    assert waveforms._common_stretch(Stream([_tr("HHZ", o, 6000), _tr("HHN", o, 6000)]), o) is None


def _cached_cli(tmp_path, monkeypatch):
    """Cached stages 1-3 and a captured orchestrator, as in test_cli_forwards_orientation_options."""
    catalog.write(catalog.from_usgs_csv(io.StringIO(COMCAT)), str(tmp_path / "event_catalog.csv"))
    (tmp_path / "network_picks.csv").touch()
    st_dir = tmp_path / "station_table"
    st_dir.mkdir()
    for name in run.STATION_FILES:
        (st_dir / name).touch()
    pd.DataFrame(dict(verdict=["used"], band=["HH"], location=[""])).to_csv(st_dir / "sensors.csv", index=False)
    called = []
    orchestrator = types.ModuleType("pocketquake.orchestrate")
    orchestrator.orchestrate = lambda *args, **kwargs: called.append((args, kwargs)) or {"status": "captured"}
    monkeypatch.setitem(sys.modules, "pocketquake.orchestrate", orchestrator)
    return called


def test_onset_events_are_the_larger_events_largest_first(tmp_path):
    catalog.write(catalog.from_usgs_csv(io.StringIO(COMCAT)), str(tmp_path / "c.csv"))
    cat = pd.read_csv(tmp_path / "c.csv")
    assert run.onset_events(cat, "2.0") == ("20151019232117", "20151019232111")
    assert run.onset_events(cat, 3.0) == ("20151019232117",)
    assert run.onset_events(cat, "none") == () and run.onset_events(cat, 5.0) == ()


@pytest.mark.parametrize("extra, expected", [
    ([], {frozenset({"20151019232117", "20151019232111"}): dict(pre=0.05, post=0.05, bandpass=(1, 40))}),
    (["--onset-mag", "3"], {frozenset({"20151019232117"}): dict(pre=0.05, post=0.05, bandpass=(1, 40))}),
    (["--onset-mag", "none"], {}),
])
def test_cluster_carries_the_onset_window(tmp_path, monkeypatch, extra, expected):
    """The generated cluster arguments build a pipeline configuration with the onset window."""
    from pipeline.clusters._base import global_cluster, layered_velmodel
    called = _cached_cli(tmp_path, monkeypatch)
    run.main(["--name", "offline", "--region", "Offline", "--work", str(tmp_path),
              "--catalog-csv", "unused.csv", "--velmodel", "con", "--apply", *extra])
    args = called[0][1]["spec_extra"]["factory_args"]
    kw = eval(f"dict({args})", {"layered_velmodel": layered_velmodel})
    cfg = global_cluster("offline", "Offline", str(tmp_path), (37.79, -121.96), (37.7, 37.9, -122.0, -121.9), **kw)
    assert cfg.xcorr_pair_overrides == expected
    assert cfg.mainshock_event_id == "20151019232117"


def test_ncsn_model_pair_follows_the_multiple_model_file(tmp_path):
    """--velmodel ncsn:gab: NCSN's primary + alternate model, delays, P-to-S ratios and alternate list
    reach HYPOINVERSE as NCSN gives them; HypoDD gets their slowness average in <= 20 layers."""
    from pipeline.clusters._base import layered_velmodel
    m = tmp_path / "models"
    m.mkdir()
    crt = lambda z, v: "X TEST 4 0.1\n" + "".join(f"{x:5.2f}" for x in z) + "\n" + "".join(f"{x:5.2f}" for x in v) + "\n"
    (m / "gab.crt").write_text(crt((0.0, 4.0, 22.5, 25.5), (3.73, 6.07, 6.47, 7.95)))
    (m / "dia.crt").write_text(crt((0.0, 1.5, 5.8, 14.4, 29.0), (2.45, 4.62, 5.80, 6.05, 7.95)))
    (m / "gab5.del").write_text("PGH   NC  0.07 D\n")
    (m / "alt.alt").write_text("BAV   NC  A\n")
    (tmp_path / "multmod.hyp").write_text(
        "* test\ndel -1 './models/alt.alt'\ncrt 5 './models/gab.crt' /pacific side\n"
        "del 5 './models/gab5.del'\ncrt 6 './models/dia.crt'\ndel 6 './models/gab5.del'\n"
        "psm 5 1.744 0\npsm 6 1.744 0\nalt 5 6\n")
    extra = tmp_path / "extra.txt"
    extra.write_text("BPV   NC  A\n")
    text, vpvs, note = run._ncsn_velocity_model_arg("ncsn_gab", "gab", str(tmp_path / "multmod.hyp"),
                                                    str(extra), str(tmp_path))
    vm = eval(text, {"layered_velmodel": layered_velmodel})
    assert vpvs == 1.744 and len(vm.p_rows) <= 20
    assert vm.hyp_lines == ("DEL -1 '{dir}/ncsn_alternate_stations.del'", "CRT 1 '{dir}/gab.crt'",
                            "DEL 1 '{dir}/gab5.del'", "CRT 2 '{dir}/dia.crt'", "DEL 2 '{dir}/gab5.del'",
                            "PSM 1 1.744 0.0", "PSM 2 1.744 0.0", "ALT 1 2")
    assert all(os.path.exists(f) for f in vm.hyp_files) and len(vm.takeoff_rows) == 2
    assert (tmp_path / "ncsn_alternate_stations.del").read_text() == "BAV   NC  A\nBPV   NC  A\n"
