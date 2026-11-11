"""Run the reduction window on the synthetic data set and summarise the result."""

import os
import hashlib
import numpy as np

import synthetic

BLOCK = 30  # pixels per side of the blocks whose sums are compared


def run_reduction(directory):
    """Run ReductionWindow (hidden, dialogs answered) on a fresh synthetic data set in `directory`."""

    synthetic.generate(directory)

    from hops import application_windows
    from hops.application_log import HOPSLog

    window_class = application_windows.HOPSWindow
    saved = {name: getattr(window_class, name) for name in ['show', 'showinfo', 'askyesno']}
    window_class.show = lambda self: self.reposition()
    window_class.showinfo = lambda self, *args, **kwargs: None
    window_class.askyesno = lambda self, *args, **kwargs: True
    try:
        from hops.application_2_reduction import ReductionWindow
        log = HOPSLog()
        synthetic.configure(log, directory)
        window = ReductionWindow(log)
        window.run(f_after=window.run_reduction)
    finally:
        for name, function in saved.items():
            setattr(window_class, name, function)

    return summarise(log, directory)


def summarise(log, directory):
    """Per reduced frame: the statistics, the sums of 30 x 30 pixel blocks and a checksum."""

    import hops.pylightcurve41 as plc
    from astropy.io import fits as pf

    all_frames = plc.open_dict(os.path.join(directory, log.all_frames))
    frames = {}
    for name in sorted(all_frames):
        frame = all_frames[name]
        data = pf.getdata(os.path.join(directory, log.reduction_directory, name), 1).astype(np.int64)
        ny, nx = data.shape[0] // BLOCK * BLOCK, data.shape[1] // BLOCK * BLOCK
        blocks = data[:ny, :nx].reshape(ny // BLOCK, BLOCK, nx // BLOCK, BLOCK).sum(axis=(1, 3))
        frames[name] = {
            'mean': float(frame[log.mean_key]),
            'std': float(frame[log.std_key]),
            'psf': float(frame[log.psf_key]),
            'jd': float(frame[log.time_key]),
            'airmass': float(frame[log.airmass_key]),
            'skip': bool(frame[log.skip_key]),
            'shape': list(data.shape),
            'block_sums': blocks.tolist(),
            'sha256': hashlib.sha256(np.ascontiguousarray(data).tobytes()).hexdigest(),
        }
    return frames
