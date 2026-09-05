"""

.. module:: test_repair_conversion
   :synopsis: tests the crosslinker-conversion figure reported after postcure repair

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import logging
logger=logging.getLogger(__name__)
from htpolynet.repair.cyanate_cap import _completion_stats

class TestCompletionStats(unittest.TestCase):
    def test_partial_cure(self):
        # a partly-cured BADCy box: 240 triazines, 58 dismantled by repair
        st=_completion_stats('TAZ',240,58)
        self.assertEqual(st['n_complete'],182)
        self.assertEqual(st['n_dismantled'],58)
        self.assertAlmostEqual(st['crosslinker_conversion'],182/240)

    def test_full_cure(self):
        st=_completion_stats('TAZ',240,0)
        self.assertEqual(st['crosslinker_conversion'],1.0)

    def test_no_crosslinkers_does_not_divide_by_zero(self):
        st=_completion_stats('TAZ',0,0)
        self.assertEqual(st['crosslinker_conversion'],0.0)

    def test_reports_the_residue_it_was_given(self):
        st=_completion_stats('XYZ',10,3)
        self.assertEqual(st['residue'],'XYZ')


class TestPreRepairBondHistogram(unittest.TestCase):
    """The histogram is the only record of how many bonds each crosslinker
    carried before repair; repair rewrites the topology and the final
    structure does not say which cap came from which ring."""

    def test_absent_when_not_supplied(self):
        st=_completion_stats('TAZ',240,58)
        self.assertNotIn('prerepair_bond_counts',st)

    def test_reported_sorted_and_int_keyed(self):
        st=_completion_stats('TAZ',240,58,log=False,
                             bond_histogram={3:182,0:5,2:40,1:13})
        self.assertEqual(list(st['prerepair_bond_counts'].keys()),[0,1,2,3])
        self.assertEqual(st['prerepair_bond_counts'],{0:5,1:13,2:40,3:182})

    def test_top_bin_is_the_surviving_population(self):
        st=_completion_stats('TAZ',240,58,log=False,
                             bond_histogram={0:5,1:13,2:40,3:182})
        self.assertEqual(st['prerepair_bond_counts'][3],st['n_complete'])

    def test_zero_bins_are_kept(self):
        # a fully cured box still reports the empty bins, so the shape of
        # repair-summary.yaml does not depend on the box
        st=_completion_stats('TAZ',240,0,log=False,
                             bond_histogram={0:0,1:0,2:0,3:240})
        self.assertEqual(st['prerepair_bond_counts'],{0:0,1:0,2:0,3:240})

    def test_disagreeing_total_warns(self):
        with self.assertLogs('htpolynet.repair.cyanate_cap',level='WARNING'):
            _completion_stats('TAZ',240,58,log=False,
                              bond_histogram={0:5,1:13,2:40,3:181})

    def test_disagreeing_top_bin_warns(self):
        # right number of residues, wrong split: 183 in the top bin against
        # 182 complete means one of the two counts is wrong
        with self.assertLogs('htpolynet.repair.cyanate_cap',level='WARNING'):
            _completion_stats('TAZ',240,58,log=False,
                              bond_histogram={0:4,1:13,2:40,3:183})

    def test_consistent_histogram_is_silent(self):
        logger_name='htpolynet.repair.cyanate_cap'
        with self.assertLogs(logger_name,level='INFO') as cm:
            _completion_stats('TAZ',240,58,
                              bond_histogram={0:5,1:13,2:40,3:182})
        self.assertFalse([r for r in cm.records if r.levelname=='WARNING'])

    def test_empty_box_does_not_warn(self):
        logger_name='htpolynet.repair.cyanate_cap'
        with self.assertLogs(logger_name,level='INFO') as cm:
            logging.getLogger(logger_name).info('probe')
            _completion_stats('TAZ',0,0,log=False,
                              bond_histogram={0:0,1:0,2:0,3:0})
        self.assertFalse([r for r in cm.records if r.levelname=='WARNING'])


class TestRunRepairWiring(unittest.TestCase):
    """The (total, stats) contract between run_repair and the runtime is only
    otherwise exercised by a full build, so pin it here."""

    def test_no_specs(self):
        from htpolynet.repair import run_repair
        total, stats = run_repair(None, None, [], None)
        self.assertEqual(total, 0)
        self.assertEqual(stats, [])

    def test_unknown_type_is_skipped_not_fatal(self):
        from htpolynet.repair import run_repair
        with self.assertLogs('htpolynet.repair', level='WARNING'):
            total, stats = run_repair(None, None, [{'type': 'no_such_driver'}], None)
        self.assertEqual(total, 0)
        self.assertEqual(stats, [])

    def test_driver_result_becomes_total_and_stats(self):
        import htpolynet.repair.cyanate_cap as cc
        from htpolynet.repair import run_repair
        recorded = _completion_stats('TAZ', 240, 58, log=False)
        original = cc.triazine_to_cyanate_cap
        cc.triazine_to_cyanate_cap = lambda *a, **k: recorded
        try:
            total, stats = run_repair(None, None, [{'type': 'triazine_to_cyanate_cap'}], None)
        finally:
            cc.triazine_to_cyanate_cap = original
        self.assertEqual(total, 58)
        self.assertEqual(stats, [recorded])
