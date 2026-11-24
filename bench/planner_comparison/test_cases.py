"""验分器的手算用例、非法序列以及冻结见证。"""

from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from cases import CASES, HERE, check, generate, is_complete


class CasesTest(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "tiny.txt"

    def fixture(self, data):
        self.path.write_text(data, encoding="ascii")
        return self.path

    def test_roster_pair_counts_and_excess(self):
        p = self.fixture("ROSTER 2 1\nS 0 1 -1\nS 0 4 7\nS 0 3 0\n")
        # 三对同时段冲突，加超上限两班。
        self.assertEqual(check(p, [0, 0, 0], []), (5, 8))
        self.assertEqual(check(p, [0, 0, 1], []), (2, 5))
        with self.assertRaises(ValueError):
            check(p, [1, 0, 0], [])
        self.assertFalse(is_complete(p, [0, 0, None], []))

    def test_timetable_conflicting_resources_add_independently(self):
        p = self.fixture("TIMETABLE 2 2\nT 0 0 1\nT 0 0 0\nT 0 1 0\n")
        # 三对教师冲突，一对班级冲突，一对教室冲突。
        self.assertEqual(check(p, [0, 0, 0, 0, 0, 1], []), (5, 1))
        with self.assertRaises(ValueError):
            check(p, [0, 0, 0, 0, 2, 0], [])

    def test_routes_close_at_depot_and_reject_duplicate(self):
        p = self.fixture("DELIVERY 2 4 1 2\nD 2 1 3\nD 3 3 2\n")
        self.assertEqual(check(p, [], [[1, 0], []]), (1, 8))
        self.assertEqual(check(p, [], [[0], [1]]), (0, 10))
        for routes in [[[0, 0], [1]], [[0], []], [[-1], [1]]]:
            self.assertFalse(is_complete(p, [], routes))
            with self.assertRaises(ValueError):
                check(p, [], routes)

    def test_factory_precedence_overlap_and_cycle(self):
        p = self.fixture("FACTORY 2\nJ 0 2 -1\nJ 1 1 0\nJ 1 2 -1\nJ 0 1 2\n")
        self.assertEqual(check(p, [], [[0, 3], [2, 1]]), (0, 3))
        self.assertEqual(check(p, [], [[3, 0], [2, 1]]), (0, 6))
        with self.assertRaisesRegex(ValueError, "循环"):
            check(p, [], [[3, 0], [1, 2]])
        with self.assertRaisesRegex(ValueError, "错误机器"):
            check(p, [], [[0, 1], [2, 3]])

    def test_frozen_witnesses(self):
        for kind, soft in zip(CASES, (719, 427, 4602, 289)):
            with self.subTest(kind=kind):
                data, values, routes = generate(kind)
                p = HERE / "instances" / (kind + ".txt")
                self.assertEqual(data, p.read_bytes())
                self.assertTrue(is_complete(p, values, routes))
                self.assertEqual(check(p, values, routes), (0, soft))


if __name__ == "__main__":
    unittest.main()
