import numpy as np
from neuroanalysis.analyzers.baseline import BaselineAnalyzer
from neuroanalysis.analyzers.stim_pulse import PatchClampStimPulseAnalyzer


class MPBaselineAnalyzer(BaselineAnalyzer):

    _settle_time = 100e-3

    @property
    def baseline_regions(self):
        """Return a list of start,stop pairs indicating regions during the recording that are expected to be quiescent
        due to absence of pulses.
        """
        if self._baseline_regions is None:
            pri = self.sync_rec.recordings[0]['primary']
            mask = np.zeros(len(pri), dtype=bool)
            dt = pri.dt
            settle_size = int(self.settle_time / dt)
            for rec in self.sync_rec.recordings:
                pa = PatchClampStimPulseAnalyzer.get(rec)
                try:
                    pulses = pa.pulses()
                except Exception:
                    print("Ignore recording baseline regions:", rec)
                    sys.excepthook(*sys.exc_info())
                    continue
                for pulse in pulses:
                    start = pri.index_at(pulse[0])
                    stop = pri.index_at(pulse[1])
                    mask[start:stop + settle_size] = True

            starts = list(np.argwhere(~mask[1:] & mask[:-1])[:,0])
            stops = list(np.argwhere(mask[1:] & ~mask[:-1])[:,0])
            if starts[0] > stops[0]:
                starts.insert(0, 0)
            if stops[-1] < starts[-1]:
                stops.append(len(mask))
            baseline_inds = [r for r in zip(starts, stops) if r[1] > r[0]]
            self._baseline_regions = [(pri.time_at(i0), pri.time_at(i1)) for i0, i1 in baseline_inds]

        return self._baseline_regions
