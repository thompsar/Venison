"""Raw Signal tab: slice averaging, phasing, and cutoff (builds the kernel)."""

from __future__ import annotations

import numpy as np
import panel as pn

from venison.core import analysis as A
from venison.ui.common import (COLORS, TIME_LABEL, Tab, button, card, grid, hint, line, plot,
                               section, sidebar_widget, vline)


class RawSignalTab(Tab):
    title = "Raw Signal"
    icon = "signal"
    step = "kernel"

    def build(self):
        spec = self.spec
        signal = A.complex_signal(spec)
        if not hasattr(spec, "phase"):
            spec.phase = A.auto_phase(signal, spec.raw_time, spec.cutoff)
        real, imag = A.phased(signal, spec.phase)

        kept = np.where(spec.slice_mask)[0]
        p_slice = plot(TIME_LABEL, height=260)
        self._slice_src = [line(p_slice, spec.raw_time, spec.real2D[kept[-1]]),
                           line(p_slice, spec.raw_time, spec.imaginary2D[kept[-1]], COLORS["imag"])]
        p_avg = plot(TIME_LABEL, height=260)
        self._avg_src = [line(p_avg, spec.raw_time, signal.real),
                         line(p_avg, spec.raw_time, signal.imag, COLORS["imag"])]
        p_phase = plot(TIME_LABEL, "Normalized echo", height=320)
        self._phase_src = [line(p_phase, spec.raw_time, real),
                           line(p_phase, spec.raw_time, imag, COLORS["imag"])]
        self._cutoff_line = vline(p_phase, spec.cutoff, color="#C0392B", dash="solid")

        self._slice_card = card(self._slice_title(kept[-1]), p_slice)
        main = [
            grid(self._slice_card, card(self._avg_title(), p_avg)),
            card("Phased signal (real, imaginary) and cutoff", p_phase),
        ]
        if self.session.is_done(self.step):
            done = "Kernel built for this phase and cutoff."
        else:
            done = "Adjust phase and cutoff, then build the kernel to continue."

        # Slices
        self._slice_pick = sidebar_widget(pn.widgets.DiscreteSlider(
            name="Slice", options=[int(i) + 1 for i in kept], value=int(kept[-1]) + 1))
        self._slice_pick.param.watch(self._on_slice, "value")
        delete = button("Delete slice")
        delete.on_click(self._on_delete)
        restore = button("Restore all slices")
        restore.on_click(self._on_restore)
        n_total = spec.real2D.shape[0]
        slices = [section("Slices"),
                  hint(f"{len(kept)} of {n_total} slices averaged.")]
        slices += [self._slice_pick, pn.Row(delete, restore, margin=0)] if n_total > 1 else []

        # Phase & cutoff
        self._phase = sidebar_widget(pn.widgets.FloatSlider(
            name="Phase (rad)", start=-np.pi, end=np.pi, step=0.001, value=float(spec.phase)))
        self._phase.param.watch(self._on_phase, "value")
        auto = button("Auto phase")
        auto.on_click(self._on_auto)
        dt = float(spec.raw_time[1] - spec.raw_time[0])
        self._cutoff = sidebar_widget(pn.widgets.FloatSlider(
            name="Cutoff (µs)", start=float(spec.raw_time[len(spec.raw_time) // 4]),
            end=float(spec.raw_time[-1]), step=dt, value=float(spec.cutoff), format="0.000"))
        self._cutoff.param.watch(self._on_cutoff, "value")
        apply = button("Build kernel", "primary")
        apply.on_click(self._on_apply)
        side = slices + [section("Phase"), self._phase, auto,
                         section("Cutoff"), self._cutoff,
                         hint("Data after the cutoff is excluded from every later step."),
                         section("Kernel"), apply, hint(done)]
        return main, side

    def _slice_title(self, idx):
        return f"Slice {idx + 1}"

    def _avg_title(self):
        return f"Average of {int(self.spec.slice_mask.sum())} slices (unphased)"

    # ── Callbacks ────────────────────────────────────────────────────────────

    def _on_slice(self, event):
        idx = int(event.new) - 1
        self._slice_src[0].data["y"] = self.spec.real2D[idx]
        self._slice_src[1].data["y"] = self.spec.imaginary2D[idx]
        self._slice_src[0].data = dict(self._slice_src[0].data)
        self._slice_src[1].data = dict(self._slice_src[1].data)
        self._slice_card[0].object = f'<div class="vn-card-title">{self._slice_title(idx)}</div>'

    def _on_delete(self, _event):
        mask = self.spec.slice_mask.copy()
        if mask.sum() <= 1:
            self.app.notify("At least one slice must stay in the average.", "warning")
            return
        mask[int(self._slice_pick.value) - 1] = False
        A.set_slice_mask(self.spec, mask)
        self.invalidate()
        self.refresh()

    def _on_restore(self, _event):
        A.set_slice_mask(self.spec, np.ones_like(self.spec.slice_mask))
        self.invalidate()
        self.refresh()

    def _update_phased(self):
        signal = A.complex_signal(self.spec)
        real, imag = A.phased(signal, self._phase.value)
        self._phase_src[0].data = {"x": self.spec.raw_time, "y": real}
        self._phase_src[1].data = {"x": self.spec.raw_time, "y": imag}

    def _on_phase(self, event):
        self.spec.phase = float(event.new)
        self._update_phased()
        self.invalidate()

    def _on_cutoff(self, event):
        self._cutoff_line.location = float(event.new)
        self.invalidate()

    def _on_auto(self, _event):
        signal = A.complex_signal(self.spec)
        self._phase.value = A.auto_phase(signal, self.spec.raw_time, self._cutoff.value)

    async def _on_apply(self, _event):
        spec = self.spec
        cutoff = spec.raw_time[A.nearest_index(spec.raw_time, self._cutoff.value)]
        phase = self._phase.value
        rebuild = getattr(spec, "tmax", None) != cutoff or not hasattr(spec, "full_kernel")
        def then(_):
            self.complete(phase=phase, cutoff=float(cutoff),
                          slices_kept=int(spec.slice_mask.sum()))
            self.app.notify("Kernel ready. Continue on the Background tab.", "success")
            self.refresh()

        await self.app.run_task(
            "Building kernel" if rebuild else "Applying phase",
            lambda: A.apply_phase_and_cutoff(spec, phase, cutoff, rebuild_kernel=rebuild), then)
