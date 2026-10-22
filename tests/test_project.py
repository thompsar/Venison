import numpy as np
import pytest

from venison.core import analysis as A
from venison.core.project import Project, create_project, list_projects
from venison.core.session import SpectrumSession

@pytest.fixture
def project(tmp_path, example_bytes):
    proj = create_project(tmp_path, "Demo")
    proj.add_spectrum("example", *example_bytes, sample="WT", condition="apo")
    return proj


def test_create_and_list(tmp_path, project):
    (summary,) = list_projects(tmp_path)
    assert summary.name == "Demo" and summary.n_spectra == 1 and summary.n_exported == 0
    with pytest.raises(FileExistsError):
        create_project(tmp_path, "Demo")
    with pytest.raises(ValueError):
        create_project(tmp_path, "../escape")


def test_labels_and_hand_dropped_files(project, example_bytes):
    (row,) = project.spectra()
    assert (row["sample"], row["condition"]) == ("WT", "apo")
    project.update_labels("example", notes="first run")
    assert project.spectra()[0]["notes"] == "first run"
    raw = project.path / "raw_data"
    (raw / "copied.DSC").write_bytes(example_bytes[0])
    (raw / "copied.DTA").write_bytes(example_bytes[1])
    assert project.spectrum_names() == ["example", "copied"]
    with pytest.raises(FileExistsError):
        project.add_spectrum("example", b"", b"")


def test_session_round_trip(project):
    s = SpectrumSession(project, "example")
    signal = A.complex_signal(s.spec)
    cutoff = s.spec.raw_time[-1]
    A.apply_phase_and_cutoff(s.spec, A.auto_phase(signal, s.spec.raw_time, cutoff), cutoff)
    s.complete("kernel")
    A.fit_zero_point(s.spec, A.symmetry_point(s.spec.raw_time))
    A.fit_background(s.spec, s.spec.background_start, "Fractal")
    s.complete("background")

    reopened = SpectrumSession(Project(project.path), "example")
    assert reopened.done == {"kernel", "background"}
    assert reopened.spec.bgmodel is A.BG_MODELS["Fractal"]
    assert isinstance(reopened.spec.cutoff, float)
    np.testing.assert_allclose(reopened.spec.waveform, s.spec.waveform)
    assert project.spectra()[0]["steps"] == ["kernel", "background"]

    reopened.complete("kernel")  # redoing an early step invalidates later ones
    assert reopened.done == {"kernel"}
