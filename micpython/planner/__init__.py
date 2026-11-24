"""Python 3.8+ 通用规划建模接口。表达式由 Rust 求值。"""
from operator import index
from threading import RLock
from weakref import finalize
from . import _native


def _int(value):
    value = index(value)
    if not -(1 << 63) <= value < (1 << 63):
        raise ValueError("整数超出 i64 范围")
    return value


def _line(*tokens):
    return " ".join(map(str, tokens))


class Expr:
    """属于一个模型的整数或布尔表达式。"""
    def __init__(self, model, node, predicate=False, var=None):
        self._model, self._node = model, node
        self._predicate, self._var = predicate, var

    def __bool__(self):
        raise TypeError("建模表达式只能由求解器判断；请使用 model.choose 或逐条 addConstr")

    def _binary(self, op, other):
        return self._model._operation(op, self, other)

    def __add__(self, other): return self._binary("add", other)
    def __radd__(self, other): return self._model._operation("add", other, self)
    def __sub__(self, other): return self._binary("sub", other)
    def __rsub__(self, other): return self._model._operation("sub", other, self)
    def __mul__(self, other): return self._binary("mul", other)
    def __rmul__(self, other): return self._model._operation("mul", other, self)
    def __eq__(self, other): return self._binary("eq", other)
    def __ne__(self, other): return self._binary("ne", other)
    def __lt__(self, other): return self._binary("lt", other)
    def __le__(self, other): return self._binary("le", other)
    def __gt__(self, other): return self._binary("gt", other)
    def __ge__(self, other): return self._binary("ge", other)


class ListVar:
    """同一组有序列表；元素下标从 0 开始。"""
    def __init__(self, model, group, owners, n):
        self._model, self._group, self.owners, self._n = model, group, owners, n

    def precede(self, before, after):
        self._model._append(_line("P", self._group, _int(before), _int(after)))

    def load(self, owner, demands):
        return self._model._node("load", self._group, _int(owner), *map(_int, demands))

    def distance(self, owner, matrix):
        rows = [list(row) for row in matrix]
        n = self._n + 1
        if len(rows) != n or any(len(row) != n for row in rows):
            raise ValueError("距离矩阵应包含各元素及最后一个仓库")
        return self._model._node("distance", self._group, _int(owner),
                                 *(_int(v) for row in rows for v in row))

    def finish(self, item):
        return self._model._node("finish", self._group, _int(item))

    def adjacent(self, owner, before, after):
        """判断指定列表中 after 是否紧跟 before；两端均为实际元素。"""
        before, after = _int(before), _int(after)
        if not (0 <= before < self._n and 0 <= after < self._n):
            raise ValueError("相邻关系的两端必须是本组元素")
        return self._model._node("adjacent", self._group, _int(owner), before, after, predicate=True)


class Result:
    """一次求解的独立快照；关闭模型后仍可读。"""
    def __init__(self, model, data):
        self._model, self._values, self._sequences = model, [], {}
        self.penalties = {}
        for line in data.splitlines():
            tag, *args = line.split()
            if tag == "S":
                self.hard, self.soft, self.steps = map(int, args[:3])
                self.status, self.complete = args[3], args[4] == "1"
            elif tag == "T":
                self.reason = args[0]
            elif tag == "V":
                self._values = [None if x == "_" else int(x) for x in args]
            elif tag == "L":
                self._sequences[int(args[0]), int(args[1])] = tuple(map(int, args[2:]))
            elif tag == "R":
                self.penalties[bytes.fromhex(args[0]).decode("utf8")] = int(args[1])

    @property
    def feasible(self):
        return self.complete and self.hard == 0

    def value(self, var):
        if not isinstance(var, Expr) or var._model is not self._model or var._var is None:
            raise ValueError("需要同一模型的分配变量")
        return self._values[var._var]

    def sequence(self, group, owner):
        if not isinstance(group, ListVar) or group._model is not self._model:
            raise ValueError("需要同一模型的列表变量")
        return self._sequences[group._group, _int(owner)]


class Model:
    """声明变量与规则，调用 optimize 返回结果。"""
    def __init__(self, name=""):
        self.name = name
        self._lines, self._vars, self._groups, self._nodes = [], 0, 0, 0
        self._names = set()
        self._objective = None
        self._closed, self._handle = False, 0
        self._failure = None
        self._lock = RLock()

    def _check(self):
        if self._closed:
            raise RuntimeError("模型已关闭")
        if self._failure is not None:
            raise ValueError(self._failure)

    def __enter__(self):
        self._check()
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        with self._lock:
            self._closed = True
            self._lines.clear()
            self._names.clear()
            if self._handle:
                self._cleanup()
                self._handle = 0

    def _append(self, line):
        self._check()
        self._lines.append(line)

    def _node(self, op, *args, predicate=False, var=None):
        self._append(_line("E", op, *args))
        expr = Expr(self, self._nodes, predicate, var)
        self._nodes += 1
        return expr

    def constant(self, value):
        return self._node("const", _int(value))

    def _expr(self, arg):
        if isinstance(arg, Expr):
            if arg._model is not self:
                raise ValueError("表达式来自另一个模型")
            return arg
        return self.constant(arg)

    def _operation(self, op, a, b):
        return self._node(op, self._expr(a)._node, self._expr(b)._node,
                          predicate=op in ("eq", "ne", "lt", "le", "gt", "ge"))

    def addVar(self, domain):
        self._append(_line("V", *map(_int, domain)))
        var = self._vars
        self._vars += 1
        return self._node("var", var, var=var)

    def addListVar(self, owners, eligible, durations=None):
        owners = _int(owners)
        eligible = [list(map(_int, row)) for row in eligible]
        durations = [0] * len(eligible) if durations is None else list(map(_int, durations))
        if len(durations) != len(eligible):
            raise ValueError("每个元素需要一个时长")
        self._append(_line("L", owners, len(eligible)))
        for d, row in zip(durations, eligible):
            self._append(_line("I", d, *row))
        group = ListVar(self, self._groups, owners, len(eligible))
        self._groups += 1
        return group

    def choose(self, condition, yes, no):
        condition = self._expr(condition)
        if not condition._predicate:
            raise TypeError("choose 需要布尔条件")
        return self._node("if", condition._node, self._expr(yes)._node, self._expr(no)._node)

    def sum(self, expressions):
        return self._node("sum", *(self._expr(x)._node for x in expressions))

    def count(self, conditions):
        return self.sum(self.choose(c, 1, 0) for c in conditions)

    def groupCount(self, keys, key):
        """统计分配给指定键的记录数。"""
        key = self._expr(key)
        return self.count(self._expr(k) == key for k in keys)

    def groupSum(self, keys, values, key):
        """累加分配给指定键的记录值。"""
        keys, values = list(keys), list(values)
        if len(keys) != len(values):
            raise ValueError("分组键与数值数量不同")
        key = self._expr(key)
        return self.sum(self.choose(self._expr(k) == key, v, 0) for k, v in zip(keys, values))

    def max(self, expressions):
        return self._node("max", *(self._expr(x)._node for x in expressions))

    def lookup(self, index_expr, table):
        return self._node("lookup", self._expr(index_expr)._node, *map(_int, table))

    def _rule(self, expr, name, weight, hard):
        if not isinstance(name, str) or not name or name in self._names:
            raise ValueError("规则名必须非空且唯一")
        weight = _int(weight)
        self._append(_line("R", name.encode("utf8").hex(), "H" if hard else "S", weight, expr._node))
        self._names.add(name)

    def addConstr(self, condition, name, weight=1, hard=True):
        self._rule(self.choose(condition, 0, 1), name, weight, hard)

    def addPenalty(self, violation, name, weight=1, hard=False):
        self._rule(self._expr(violation), name, weight, hard)

    def setObjective(self, expr):
        self._check()
        self._objective = self._expr(expr)._node

    def optimize(self, seed=42, steps=50000, time_ms=-1, threads=1, algo="tabu", construct="ff"):
        with self._lock:
            self._check()
            seed = index(seed)
            if not 0 <= seed < (1 << 64):
                raise ValueError("seed 须在 u64 范围内")
            if algo not in ("tabu", "sa", "la") or construct not in ("ff", "bf"):
                raise ValueError("未知搜索算法或构造方式")
            extra = [] if self._objective is None else [_line("O", self._objective)]
            extra.append(_line("C", seed, _int(steps), _int(time_ms), _int(threads), algo, construct))
            request = "\n".join(self._lines + extra).encode("ascii")
            lib = _native.library()
            if not self._handle:
                self._handle = lib.planner_create()
                if not self._handle:
                    raise RuntimeError("建模句柄创建失败")
                self._cleanup = finalize(self, lib.planner_destroy, self._handle)
            try:
                response = _native.run(lib, self._handle, request)
            except ValueError as error:
                self._failure = str(error)
                raise
            return Result(self, response)


__all__ = ["Model", "Expr", "ListVar", "Result"]
