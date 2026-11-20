"""将固定版本的 PyBoolNet 原规则转成定点 MILP，并用原规则穷举对账。"""
import argparse
import ast
import hashlib
import itertools
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path


SOURCE = "3394310d07e3e25810a828d8da4f6d7241115113"
DIRECTORY = "pyboolnet/repository"
OPTIONS = ["--threads", "1", "--node-limit", "2000000", "--time-limit-secs", "60",
           "--early-stop", "off"]
TOKEN = re.compile(r"[A-Za-z_][A-Za-z_0-9]*|[01!&|()]")


def read_rules(raw):
    rules = {}
    for line in raw.decode().splitlines():
        line = line.partition("#")[0].strip()
        if not line or line.replace(" ", "") == "targets,factors":
            continue
        name, expr = map(str.strip, line.split(","))
        if name in rules:
            raise ValueError(f"重复定义 {name}")
        tokens = TOKEN.findall(expr)
        if "".join(tokens) != "".join(expr.split()):
            raise ValueError(f"未知符号 {expr}")
        rules[name] = tokens
    return rules


def equations(rules):
    # 只允许上游固定快照中的布尔符号；先映射变量名，避免与 Python 保留字冲突。
    names = list(rules)
    symbols = dict(zip(names, (f"x{i}" for i in range(len(names)))))
    translation = {**symbols, "&": "and", "|": "or", "!": "not",
                   "0": "0", "1": "1", "(": "(", ")": ")"}
    return [" ".join(translation[token] for token in rhs) for rhs in rules.values()]


def fixed_points(expressions):
    # 直接按 Python 布尔运算执行原始规则，此处不读 MILP 行或辅助变量。
    variables = ",".join(f"x{i}" for i in range(len(expressions)))
    rhs = ",".join(f"({expression})" for expression in expressions)
    update = eval(f"lambda {variables}: ({rhs},)", {"__builtins__": {}})
    return {state for state in itertools.product((0, 1), repeat=len(expressions))
            if state == update(*state)}


class Model:
    def __init__(self, name, expressions):
        self.name = name
        self.columns = len(expressions)
        self.rows = []
        for i, expression in enumerate(expressions):
            self.define((i, True), ast.parse(expression, mode="eval").body)

    def row(self, terms, op, rhs):
        # (列, 正文字) 表示 x，(列, 负文字) 表示 1-x；常量也可以直接移到右端。
        coefficients = Counter()
        for weight, (column, positive) in terms:
            if column is None:
                rhs -= weight * int(positive)
            elif positive:
                coefficients[column] += weight
            else:
                coefficients[column] -= weight
                rhs -= weight
        row = sorted((col, value) for col, value in coefficients.items() if value)
        if not row and (rhs == 0 if op == "=" else rhs >= 0):
            return
        self.rows.append((row, op, rhs))

    def literal(self, node):
        if isinstance(node, ast.Name):
            return int(node.id[1:]), True
        if isinstance(node, ast.Constant):
            return None, bool(node.value)
        if isinstance(node, ast.UnaryOp):
            column, positive = self.literal(node.operand)
            return column, not positive
        column = self.columns
        self.columns += 1
        self.define((column, True), node)
        return column, True

    def define(self, output, node):
        if isinstance(node, ast.UnaryOp):
            self.define((output[0], not output[1]), node.operand)
        elif isinstance(node, ast.BoolOp):
            inputs = [self.literal(child) for child in node.values]
            width = len(inputs)
            # 整数域的两条门约束。与门满足 k*y <= sum <= k-1+y，
            # 或门满足 y <= sum <= k*y。根门直接使用原变量，不添加复制层。
            if isinstance(node.op, ast.And):
                lower, upper, offset = width, 1, width - 1
            else:
                lower, upper, offset = 1, width, 0
            self.row([(lower, output)] + [(-1, lit) for lit in inputs], "<=", 0)
            self.row([(-upper, output)] + [(1, lit) for lit in inputs], "<=", offset)
        else:
            self.row([(1, output), (-1, self.literal(node))], "=", 0)

    def save(self, target):
        lines = [f"name {self.name}", "sense min", "vars"]
        lines += [f"x{i} binary" for i in range(self.columns)]
        lines += ["obj", "1 x0", "cons"]
        for row, op, rhs in self.rows:
            lhs = " ".join(f"{value} x{column}" for column, value in row)
            lines.append(f"{lhs} {op} {rhs}")
        target.write_text("\n".join(lines + ["end", ""]))

    def verify(self, event, points, n):
        values = event["values"]
        if values is None:
            return event["objective"] is None
        if len(values) != self.columns or any(abs(x - round(x)) > 1e-6 or
                                               round(x) not in (0, 1) for x in values):
            return False
        state = tuple(round(x) for x in values)
        if state[:n] not in points or abs(event["objective"] - state[0]) > 1e-6:
            return False
        for row, op, rhs in self.rows:
            lhs = sum(state[i] * coeff for i, coeff in row)
            if (lhs != rhs if op == "=" else lhs > rhs):
                return False
        return True


def verify(record, model, points, n):
    events = record.get("incumbent_events", [])
    target = min((x[0] for x in points), default=None)
    if points:
        correct = (record["status"] == "optimal" and record["values"] is not None and
                   record["objective"] is not None and abs(record["objective"] - target) <= 1e-6 and
                   record["bound"] is not None and abs(record["bound"] - target) <= 1e-6)
    else:
        correct = record["status"] == "infeasible" and record["values"] is None
    verified = (correct and not record.get("callback_errors") and model.verify(record, points, n)
                and all(e["values"] is not None and model.verify(e, points, n) for e in events))
    def bits(event):
        return "".join(str(round(v)) for v in event["values"]) if event["values"] is not None else None

    return {**{k: record[k] for k in ("status", "objective", "bound", "nodes", "lp_iterations",
                                      "elapsed_secs")}, "callbacks": len(events),
            "final_bits": bits(record), "callback_bits": [bits(e) for e in events],
            "verified": verified}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for opt in ("repo", "models", "baseline", "current", "output", "detector"):
        parser.add_argument(f"--{opt}", type=Path, required=True)
    parser.add_argument("--gurobi", type=Path)
    args = parser.parse_args()
    args.models.mkdir(exist_ok=True, parents=True)
    if args.gurobi:
        sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
        from compare_drayage_timing import observe_gurobi
        from compare_gurobi10 import solve

    def git(*params):
        return subprocess.check_output(["git", "-C", args.repo, *params])

    # 直接从固定 Git 对象读取，也适用于工作区存在本地修改的情况。
    paths = [Path(f) for f in git("ls-tree", "-r", "--name-only", SOURCE, DIRECTORY)
             .decode().splitlines() if f.endswith(".bnet")]
    seen, skipped, results = {}, {}, []
    for path in paths:
        name = path.stem
        if name != path.parent.name:
            skipped[str(path)] = "非同名主模型"
            continue
        raw = git("show", f"{SOURCE}:{path.as_posix()}")
        rules = read_rules(raw)
        if not 2 <= len(rules) <= 16:
            skipped[name] = f"{len(rules)} 个变量超出 2–16"
            continue
        if name == "multivalued":
            skipped[name] = "上游注明采用多值水平编码"
            continue
        key = tuple(sorted((name, tuple(tokens)) for name, tokens in rules.items()))
        if key in seen:
            skipped[name] = f"与 {seen[key]} 的原规则相同"
            continue
        seen[key] = name
        expr = equations(rules)
        points = fixed_points(expr)
        model = Model(name, expr)
        dest = args.models / f"{name}.milp"
        model.save(dest)
        row = {"name": name, "source_sha256": hashlib.sha256(raw).hexdigest(),
               "file_sha256": hashlib.sha256(dest.read_bytes()).hexdigest(),
               "states": 2 ** len(rules), "state_variables": len(rules),
               "objective_variable": next(iter(rules)), "fixed_points": len(points),
               "minimum": min((x[0] for x in points), default=None),
               "milp_variables": model.columns, "milp_rows": len(model.rows),
               "two_column_rows": sum(len(terms) == 2 for terms, _, _ in model.rows)}
        row["detector"] = json.loads(subprocess.check_output([args.detector, dest]))
        runs = {}
        for variant in ("baseline", "current"):
            command = [getattr(args, variant), "timed-solve", dest, *OPTIONS]
            runs[variant] = lambda cmd=command: json.loads(subprocess.check_output(cmd, timeout=75))
        if args.gurobi:
            data = json.loads(subprocess.check_output([args.current, "export", dest]))
            runs["gurobi"] = lambda: solve(args.gurobi, data, 60, 2_000_000, observe_gurobi)
        for variant, run in runs.items():
            row[variant] = verify(run(), model, points, len(rules))
        if row["detector"]["eligible"]:
            # 单次已预热；命中样例再轮换三方的先后次序，顺序运行并逐次核验。
            names = list(runs)
            timings = {v: [] for v in names}
            for i in range(9):
                start = i % len(names)
                for variant in names[start:] + names[:start]:
                    trial = verify(runs[variant](), model, points, len(rules))
                    if not trial["verified"]:
                        raise ValueError(f"{name} 重跑校验失败 {variant}: {trial}")
                    timings[variant].append(trial["elapsed_secs"])
            row["verified_repeat_secs"] = timings
        results.append(row)
    audit = {"source_revision": SOURCE, "probe_options": OPTIONS, "source_files": len(paths),
             "skipped": skipped, "results": results}
    args.output.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n")
    variants = ("baseline", "current", "gurobi") if args.gurobi else ("baseline", "current")
    if not results or not all(r[v]["verified"] for r in results for v in variants):
        raise SystemExit("有记录未通过定点、证明状态或原约束核验，请检查结果")


if __name__ == "__main__":
    main()
