"""用固定矩阵的合成汇总验证单特性门禁，无需求解器。"""
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from gate_feature_trial import decide, main
from run_feature_trial import FEATURES, VARIANTS, expected_features
from run_root_work import CASES
from summarize_feature_trial import MIN_RSS_LIMIT_BYTES, summarize


def make_summary(feature=FEATURES[0], phase='paired', ratio=1.0, stable=0, *,
                 timed_models=3, mutate=None, errors=()):
    seconds, reps = (10, 1) if phase == 'screen' else (30, 3)
    rows = []
    for index, name in enumerate(CASES):
        for variant in VARIANTS:
            for rep in range(reps):
                optimal = index < timed_models
                if variant == 'candidate':
                    optimal |= timed_models <= index < timed_models + stable
                row = {
                    'instance': name, 'variant': variant, 'rep': rep,
                    'status': 'optimal' if optimal else 'limit_reached',
                    'objective': 1.0 if optimal else None, 'bound': 1.0 if optimal else None,
                    'values': [1.0] if optimal else None, 'elapsed_secs': 1.0 if variant == 'base' else ratio,
                    'auto_degraded': False, 'gap_terminated': False, 'peak_rss_bytes': 1024,
                    'nodes': 1, 'lp_iterations': 1,
                }
                if mutate is not None:
                    mutate(row, index)
                rows.append(row)
    source = dict(source_commit='合成提交', source_tree='合成源码树', binary_manifest_sha256='合成摘要', host='合成主机')
    builds = {variant: dict(binary=variant, binary_sha256=variant, features=expected_features(feature, variant))
              for variant in VARIANTS}
    checked = dict(passed=not errors, errors=list(errors), runs=len(rows), instances=len(CASES))
    with (patch('summarize_feature_trial.load_runs', return_value=(rows, feature, phase, seconds, reps, source, builds)),
          patch('summarize_feature_trial.inspect_results', return_value=checked)):
        return summarize(Path('合成记录.jsonl'))


class GateFeatureTrialTest(unittest.TestCase):
    def test_paired_performance_boundaries(self):
        for ratio, stable, passed in ((0.9, 0, True), (0.9001, 0, False), (1.05, 2, True),
                                      (1.0501, 2, False), (1.05, 1, False)):
            with self.subTest(ratio=ratio, stable=stable):
                result = decide(make_summary(ratio=ratio, stable=stable))
                self.assertIs(result['eligible_for_default'], passed)
                self.assertEqual(result['new_stable_proof_models'], list(CASES[3:3 + stable]))

    def test_paired_requires_three_eligible_cases(self):
        result = decide(make_summary(ratio=0.8, timed_models=2))
        self.assertTrue(result['gates']['performance'])
        self.assertFalse(result['eligible_for_default'])

    def test_paired_requires_no_loss_even_if_candidate_proves_more_models(self):
        def move_one_proof(row, index):
            if index == 0 and row['variant'] == 'candidate' and row['rep'] == 1:
                row.update(status='limit_reached')

        result = decide(make_summary(ratio=0.8, stable=2, timed_models=4, mutate=move_one_proof))
        self.assertEqual(result['new_stable_proof_models'], list(CASES[4:6]))
        self.assertEqual(result['paired_proof_losses'], 1)
        self.assertFalse(result['eligible_for_default'])

    def test_base_partial_proof_is_not_new_stable_proof(self):
        def prove_once_in_base(row, index):
            if index == 3 and row['variant'] == 'base' and row['rep'] == 1:
                row.update(status='optimal', objective=1.0, bound=1.0, values=[1.0])

        result = decide(make_summary(stable=2, mutate=prove_once_in_base))
        self.assertEqual(result['new_stable_proof_models'], list(CASES[4:5]))
        self.assertFalse(result['eligible_for_default'])

    def test_screen_rss_and_case_time_gates_depend_on_feature(self):
        def break_one_gate(row, index, gate):
            if index == 0:
                row['elapsed_secs'] = 0.03125
                if row['variant'] == 'candidate':
                    if gate == 'rss':
                        row['peak_rss_bytes'] = MIN_RSS_LIMIT_BYTES + 1
                    else:
                        row['elapsed_secs'] = 0.0625

        for feature in FEATURES:
            strict = feature in ('root_heuristic_work_limit', 'certified_row_fallback')
            for gate in ('rss', 'case_median'):
                with self.subTest(feature=feature, gate=gate):
                    result = decide(make_summary(feature, 'screen',
                                                 mutate=lambda row, index: break_one_gate(row, index, gate)))
                    self.assertEqual(result['timing']['models_over_case_ratio'],
                                     list(CASES[:1]) if gate == 'case_median' else [])
                    self.assertEqual(result['timing']['eligible_models'], 2)
                    self.assertEqual(result['rss_failures'], int(gate == 'rss'))
                    self.assertIs(result['eligible_for_paired'], not strict)
                    if strict:
                        self.assertFalse(result['gates'][gate])

            with self.subTest(feature=feature, ratio=1.25):
                result = decide(make_summary(feature, 'screen', ratio=1.25))
                self.assertTrue(result['eligible_for_paired'])

    def test_correctness_failure_blocks_promotion(self):
        result = decide(make_summary(ratio=0.8, errors=['合成正确性失败']))
        self.assertFalse(result['gates']['correctness'])
        self.assertFalse(result['eligible_for_default'])

    def test_rejects_conflicting_optimal_objectives_in_new_proof(self):
        summary = make_summary(ratio=1.0, stable=2)
        candidate = summary['cases'][CASES[3]]['pairs'][1]['candidate']
        candidate.update(objective=2.0, bound=2.0)
        with self.assertRaisesRegex(ValueError, '已证最优目标.*不一致'):
            decide(summary)

    def test_rejects_missing_or_inconsistent_evidence(self):
        for edit in (
                lambda s: s.pop('source'),
                lambda s: s.update(repetitions=True),
                lambda s: s.update(phase='screen'),
                lambda s: s.update(paired_proof_losses=1),
                lambda s: s['variants']['base'].update(proof_runs=0),
                lambda s: s['timing'].update(eligible_models=0),
                lambda s: s['timing'].update(geometric_candidate_over_base=float('inf')),
                lambda s: s['cases'][CASES[0]]['pairs'][0].update(rep=1),
                lambda s: s['cases'][CASES[0]].update(common_optimal_timing=None),
                lambda s: s['cases'][CASES[0]]['pairs'][0].update(rss_passed=False)):
            summary = make_summary()
            edit(summary)
            with self.subTest(summary_edit=edit), self.assertRaises((KeyError, ValueError)):
                decide(summary)

    def test_cli_rejects_nonfinite_and_duplicate_json(self):
        payload = json.dumps(make_summary())
        for suffix in ('"extra": NaN}', '"extra": Infinity}', '"extra": 1e400}',
                       '"correctness_passed": true}'):
            with self.subTest(suffix=suffix), TemporaryDirectory() as tmp:
                summary = Path(tmp) / 'summary.json'
                summary.write_text(payload[:-1] + ', ' + suffix, encoding='utf-8')
                with (patch('sys.argv', ['gate_feature_trial.py', str(summary)]),
                      redirect_stdout(StringIO()), redirect_stderr(StringIO())):
                    self.assertEqual(main(), 2)


if __name__ == '__main__':
    unittest.main()
