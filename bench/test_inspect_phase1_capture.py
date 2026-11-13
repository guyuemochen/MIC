import json
from pathlib import Path
import struct
import unittest

from inspect_phase1_capture import inspect


FIXTURES = Path(__file__).parent / 'large_scale_rounds/results/phase1-progress-validation'


def bits(value):
    return struct.pack('>d', value).hex()


def fixture(engine='revised'):
    return json.loads((FIXTURES / f'06-capture-{engine}.json').read_bytes())


def candidate(document):
    key = 'snapshots'
    snapshots = document[key] if key in document else document['capture'][key]
    return next(snapshot for snapshot in snapshots if snapshot['kind'] == 'first_candidate')


class InspectPhase1CaptureTest(unittest.TestCase):
    def test_real_fixture_and_same_state_pricing(self):
        for engine in ('dense', 'revised'):
            with self.subTest(engine=engine):
                doc = fixture(engine)
                # 外部坐标不能覆盖已采集的引擎模型。
                doc['external_input_not_engine_coordinates']['rows'][0]['rhs'] = bits(-10)
                report = inspect(doc)
                self.assertEqual(report['issues'], [])
                self.assertEqual(report['candidate_conclusion'], 'primal_and_augmented_pass')
                self.assertEqual(report['same_state_candidate_pricing'][0]['kind'], 'optimal')
                self.assertFalse(report['feasible_candidate_with_entering_column'])
                pricing = next(e for e in doc['capture']['events'] if e['kind'] == 'pricing_at_first_candidate')
                pricing['pricing'].update(kind='candidate', column=1, direction=bits(1.0))
                self.assertTrue(inspect(doc)['feasible_candidate_with_entering_column'])
                pricing['pricing']['evaluated_revision'] -= 1
                self.assertFalse(inspect(doc)['feasible_candidate_with_entering_column'])

    def test_current_bounds_and_rows_are_distinct(self):
        doc = fixture()
        candidate(doc)['upper'][0] = bits(999.0)
        result = candidate(inspect(doc))
        self.assertEqual(result['conclusion'], 'near_zero_but_primal_infeasible')
        self.assertEqual(result['current_structural_bounds']['max_violation'], 1.0)
        self.assertTrue(result['construction_bounds']['passed'])
        self.assertTrue(result['rows']['normalized']['passed'])

        doc = fixture()
        snapshot = candidate(doc)
        snapshot['basic_values'][0] = snapshot['x_unclipped'][0] = bits(999.0)
        result = candidate(inspect(doc))
        self.assertTrue(result['current_structural_bounds']['passed'])
        self.assertEqual(result['rows']['normalized']['max_violation'], 1.0)
        self.assertEqual(result['rows']['engine_input']['raw_max_violation'], 1.0)
        self.assertEqual(result['rows']['normalized']['worst']['source'], {'kind': 'constraint', 'index': 0})

    def test_artificial_location_is_not_basis_position(self):
        doc = fixture()
        snapshot = candidate(doc)
        matrix = doc['capture']['normalized_matrix']
        # 人工列编号不变，调整原行归属并使用负号。
        matrix['artificial_columns'] = [[1, bits(1)], [0, bits(-1)]]
        value = bits(2e-10)
        snapshot['basic_values'][1] = value
        snapshot['summary'].update(positive_sum=value, minimum=value, maximum=value, artificial_nonzero=1)
        result = candidate(inspect(doc))
        self.assertTrue(result['artificial']['near_zero'])
        self.assertEqual([row['absolute_residual'] for row in result['rows']['augmented']['rows']], [2e-10, 0.0])
        snapshot['factorization_broken'] = True
        self.assertEqual(inspect(doc)['candidate_conclusion'], 'primal_pass_but_augmented_inconsistent')

    def test_missing_truncated_and_invalid_mapping(self):
        doc = fixture()
        doc['capture']['snapshots'].remove(candidate(doc))
        doc['capture']['slot_omitted'][1] = True
        report = inspect(doc)
        self.assertEqual(report['capture']['slots']['first_candidate']['state'], 'omitted')
        self.assertEqual(report['candidate_conclusion'], 'insufficient_evidence')
        doc['capture']['slot_omitted'][1] = False
        doc['capture']['slot_attempt_iterations'][1] = None
        self.assertEqual(inspect(doc)['capture']['slots']['first_candidate']['state'], 'not_triggered')
        doc = fixture()
        doc['capture']['normalized_matrix']['artificial_columns'][0][1] = bits(2)
        self.assertIn('artificial_columns', inspect(doc)['issues'][0])
        doc = fixture()
        doc['capture']['engine_input_lp'] = None
        self.assertIn('engine_input_lp', inspect(doc)['issues'][0])
        doc = fixture()
        candidate(doc)['x_unclipped'][0] = '0.0'
        self.assertIn('16 位十六进制', inspect(doc)['issues'][0])

    def test_nonfinite_payloads_illegal_upper_and_negative_zero(self):
        for payload in ('7ff8000000000123', '7ff0000000000000', 'fff0000000000000'):
            with self.subTest(payload=payload):
                doc = fixture()
                snapshot = candidate(doc)
                snapshot['basic_values'][0] = snapshot['x_unclipped'][0] = payload
                result = candidate(inspect(doc))
                self.assertEqual(result['conclusion'], 'insufficient_evidence')
                self.assertEqual(result['nonfinite_reconstructed_columns'], [{'column': 0, 'bits': payload}])
        doc = fixture()
        candidate(doc)['at_upper'][2] = True
        result = candidate(inspect(doc))
        self.assertEqual(result['nonfinite_reconstructed_columns'][0]['column'], 2)
        doc = fixture()
        candidate(doc)['x_unclipped'][1] = bits(-0.0)
        result = candidate(inspect(doc))
        self.assertEqual(result['x_bit_mismatches'][0]['saved'], bits(-0.0))
        self.assertEqual(result['conclusion'], 'insufficient_evidence')

    def test_input_near_zero_and_unscaled_row_differences(self):
        doc = fixture()
        capture = doc['capture']
        matrix = capture['normalized_matrix']
        model = capture['engine_input_lp']
        model['rows'][0]['terms'] = [[0, bits(4)], [1, bits(4)]]
        model['rows'][0]['rhs'] = bits(4000)
        matrix['row_scale'][0] = bits(0.25)
        model['rows'].append({'terms': [[0, bits(1e-12)]], 'rhs': bits(-1e-9), 'op': 'Le'})
        matrix['constraint_row'].append(None)
        result = candidate(inspect(doc))
        self.assertTrue(result['rows']['engine_input']['scaled']['passed'])
        self.assertFalse(result['rows']['engine_input']['dropped_near_zero']['passed'])
        self.assertEqual(result['rows']['engine_input']['raw_worst']['constraint'], 2)

        # 同一行按当前约束来源变成割，不应套用输入行的缩放规则。
        model['rows'].pop(0)
        matrix['constraint_row'].pop(0)
        matrix['cut_row'] = [0]
        matrix['row_scale'][0] = bits(1.0)
        snapshot = candidate(doc)
        snapshot['basic_values'][0] = snapshot['x_unclipped'][0] = bits(999.0)
        result = candidate(inspect(doc))
        self.assertEqual(result['rows']['cuts']['worst']['source'], {'kind': 'cut', 'index': 0})
        self.assertEqual(result['rows']['cuts']['max_violation'], 1.0)

    def test_overflow_is_retained_and_cancellation_is_order_independent(self):
        doc = fixture()
        snapshot = candidate(doc)
        matrix = doc['capture']['normalized_matrix']
        # 只扰动引擎原尺度。归一矩阵正确，原乘积和 RHS 的绝对和越过双精度范围。
        doc['capture']['engine_input_lp']['rows'][0]['terms'] = [[0, bits(1e305)], [1, bits(1e305)]]
        doc['capture']['engine_input_lp']['rows'][0]['rhs'] = bits(1e308)
        matrix['row_scale'][0] = bits(1.0 / 1e305)
        matrix['row_rhs'][0] = bits(1e308 / 1e305)
        snapshot['basic_values'][0] = snapshot['x_unclipped'][0] = matrix['row_rhs'][0]
        for reverse in (False, True):
            if reverse:
                doc['capture']['engine_input_lp']['rows'][0]['terms'].reverse()
            result = candidate(inspect(doc))
            self.assertEqual(result['conclusion'], 'insufficient_evidence')
            self.assertTrue(result['rows']['engine_input']['rows'][0]['exact_recheck'])
            self.assertEqual(result['rows']['binary64_arithmetic_issues'][0]['layer'], 'engine_input')
            raw = result['rows']['engine_input']['rows'][0]['raw_violation']
            if not reverse:
                original_raw = raw
            else:
                self.assertEqual(raw, original_raw)


if __name__ == '__main__':
    unittest.main()
