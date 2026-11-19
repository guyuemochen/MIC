import unittest

from compare_drayage_timing import check_events


class IncumbentObservationTest(unittest.TestCase):
    def test_rejects_bad_event_and_still_reports_first_route(self):
        model = dict(sense='min', offset=0, objective=[[0, 5.0]],
                     variables=[dict(kind='binary', name='skip', lb=0.0, ub=1.0)], rows=[])
        row = dict(status='optimal', elapsed_secs=4.0, objective=0.0, values=[0.0],
                   incumbent_events=[
                       dict(elapsed_secs=1.0, objective=0.0, values=[-1.0]),
                       dict(elapsed_secs=2.0, objective=5.0, values=[1.0]),
                       dict(elapsed_secs=3.0, objective=0.0, values=[0.0]),
                   ])
        observed = check_events(model, row, [0])
        self.assertEqual(observed['errors'], [])
        self.assertEqual(observed['rejected'][0]['event'], 1)
        self.assertEqual(observed['first']['event'], 2)
        self.assertEqual(observed['first']['tasks_served'], 0)
        self.assertEqual(observed['first_route']['event'], 3)
        self.assertEqual(observed['first_route']['tasks_served'], 1)
        self.assertEqual(observed['final_vector_event'], 3)

    def test_final_vector_is_not_a_first_solution_event(self):
        model = dict(sense='min', offset=0, objective=[], variables=[], rows=[])
        row = dict(status='optimal', elapsed_secs=4.0, objective=0.0, values=[])
        observed = check_events(model, row)
        self.assertIsNone(observed['first'])
        self.assertTrue(observed['errors'])
