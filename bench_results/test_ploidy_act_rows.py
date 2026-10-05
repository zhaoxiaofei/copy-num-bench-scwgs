#!/usr/bin/env python3
"""Regression test for the ACT row identity in the ploidy balloon figure.

TN6 and TN7 were each sequenced at two read lengths (36 bp and 152 bp), and
the two runs differ in both the number of evaluable cells and the percentage
inside the ploidy window.  Before this test existed the ACT panel keyed its
grid cells by the bare sample id, so the last-loaded run silently won: the
figure showed one ambiguous "TN6" row whose cell count could be read as
either run, while the long table kept both and the statistics aggregated
them.  These checks lock in the fix:

  1. ``act_row_label()`` qualifies every ACT row with its read length;
  2. ``order_act_rows()`` keeps the two runs of one sample adjacent and
     ordered by read length;
  3. ``_build_payload()`` never silently drops one of two same-key runs: it
     keeps the more evaluable row and logs the collision.
"""
import importlib.util
import logging
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, 'scWGS-ploidy-performances-eval.py')

spec = importlib.util.spec_from_file_location('ploidy_eval_script', SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def _entry(dataset, tool, cap, n, pct, spot=''):
    return {'plot': 'ACT', 'dataset': dataset, 'tool': tool, 'max_cn': cap,
            'method': F'{tool}|{mod.fmt_max_cn(cap)}', 'window': 0.5,
            'n_cells': n, 'n_cells_finite': n, 'pct_within': pct,
            'mean_abs_ploidy_error': 0.1, 'expected_ploidy_mean': 3.0,
            'failed': False, 'donor': 'D', 'sampleType': 'tumor',
            'avgSpotLen': spot, 'cellLine': ''}


class _Capture(logging.Handler):
    def __init__(self):
        super().__init__()
        self.messages = []

    def emit(self, record):
        self.messages.append(record.getMessage())


def main():
    # 1. every ACT row label carries its read length.
    assert mod.act_row_label('TN6', '36') == 'TN6 · 36 bp'
    assert mod.act_row_label('TN6', '152') == 'TN6 · 152 bp'
    assert mod.act_row_label('TN1', '') == 'TN1'
    print('OK 1: ACT row labels carry avgSpotLen')

    # 2. TN6/TN7's two runs stay adjacent and ordered by read length.
    sub = pd.DataFrame({'dataset': ['TN1 · 50 bp', 'TN6 · 152 bp', 'TN6 · 36 bp',
                                    'TN7 · 152 bp', 'TN7 · 36 bp', 'mb157 · 50 bp']})
    assert mod.order_act_rows(sub) == ['TN1 · 50 bp', 'TN6 · 36 bp',
                                       'TN6 · 152 bp', 'TN7 · 36 bp',
                                       'TN7 · 152 bp', 'mb157 · 50 bp']
    print('OK 2: TN6/TN7 runs stay adjacent and read-length ordered')

    # 3a. distinguishable runs (the real TN6 case) keep their own grid cells.
    entries = pd.DataFrame([
        _entry('TN6 · 152 bp', 'ginkgo', 10.0, 173, 93.6, spot='152'),
        _entry('TN6 · 36 bp', 'ginkgo', 10.0, 1205, 99.1, spot='36'),
    ])
    payload = mod._build_payload(entries, 'ACT',
                                 ['TN6 · 36 bp', 'TN6 · 152 bp'],
                                 [('ginkgo', 10.0)])
    assert payload['count'][0, 0] == 1205
    assert payload['count'][1, 0] == 173
    print('OK 3a: TN6 36 bp and 152 bp are separate grid cells')

    # 3b. a true same-key collision is never resolved silently: the more
    #     evaluable row is kept and the collision is logged.
    handler = _Capture()
    root = logging.getLogger()
    root.addHandler(handler)
    try:
        dup = pd.DataFrame([
            _entry('TN6 · 36 bp', 'ginkgo', 10.0, 173, 50.0, spot='36'),
            _entry('TN6 · 36 bp', 'ginkgo', 10.0, 1205, 90.0, spot='36'),
        ])
        payload = mod._build_payload(dup, 'ACT', ['TN6 · 36 bp'],
                                     [('ginkgo', 10.0)])
    finally:
        root.removeHandler(handler)
    assert payload['count'][0, 0] == 1205
    assert any('two evaluations match' in m for m in handler.messages), \
        handler.messages
    print('OK 3b: same-key collision keeps the more evaluable row and warns')

    print('\nAll ACT row-identity regression checks PASSED')


if __name__ == '__main__':
    main()
