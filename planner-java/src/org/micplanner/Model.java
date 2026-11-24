package org.micplanner;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/** 单个业务模型。建模由调用方串行完成，各个模型可以并行求解。 */
public final class Model implements AutoCloseable {
    private final StringBuilder lines = new StringBuilder();
    private final Set<String> names = new HashSet<String>();
    private int vars, groups, nodes;
    private Expr objective;
    private String failure;
    private boolean closed;

    private void check() {
        if (closed) throw new IllegalStateException("模型已关闭");
        if (failure != null) throw new IllegalStateException(failure);
    }

    private void line(String line) {
        check();
        lines.append(line).append('\n');
    }

    private static void longs(StringBuilder b, long[] values) {
        for (long value : values) b.append(' ').append(value);
    }

    private static String hex(String name) {
        StringBuilder out = new StringBuilder();
        for (byte b : name.getBytes(StandardCharsets.UTF_8)) {
            out.append(Character.forDigit((b & 255) >>> 4, 16));
            out.append(Character.forDigit(b & 15, 16));
        }
        return out.toString();
    }

    private static String unhex(String text) {
        byte[] bytes = new byte[text.length() / 2];
        for (int i = 0; i < bytes.length; i++) {
            bytes[i] = (byte)Integer.parseInt(text.substring(i * 2, i * 2 + 2), 16);
        }
        return new String(bytes, StandardCharsets.UTF_8);
    }

    private Expr node(String op, boolean predicate, int var) {
        line("E " + op);
        return new Expr(this, nodes++, predicate, var);
    }

    private int id(Expr e) {
        if (e == null || e.model != this) throw new IllegalArgumentException("表达式来自另一个模型");
        return e.node;
    }

    public static final class Expr {
        private final Model model;
        private final int node, var;
        private final boolean predicate;

        private Expr(Model model, int node, boolean predicate, int var) {
            this.model = model;
            this.node = node;
            this.predicate = predicate;
            this.var = var;
        }

        private Expr binary(String op, Expr other, boolean predicate) {
            return model.node(op + " " + node + " " + model.id(other), predicate, -1);
        }

        public Expr add(Expr b) { return binary("add", b, false); }
        public Expr sub(Expr b) { return binary("sub", b, false); }
        public Expr mul(Expr b) { return binary("mul", b, false); }
        public Expr eq(Expr b) { return binary("eq", b, true); }
        public Expr ne(Expr b) { return binary("ne", b, true); }
        public Expr lt(Expr b) { return binary("lt", b, true); }
        public Expr le(Expr b) { return binary("le", b, true); }
        public Expr gt(Expr b) { return binary("gt", b, true); }
        public Expr ge(Expr b) { return binary("ge", b, true); }

        public Expr add(long b) { return add(model.constant(b)); }
        public Expr sub(long b) { return sub(model.constant(b)); }
        public Expr mul(long b) { return mul(model.constant(b)); }
        public Expr eq(long b) { return eq(model.constant(b)); }
        public Expr ne(long b) { return ne(model.constant(b)); }
        public Expr lt(long b) { return lt(model.constant(b)); }
        public Expr le(long b) { return le(model.constant(b)); }
        public Expr gt(long b) { return gt(model.constant(b)); }
        public Expr ge(long b) { return ge(model.constant(b)); }
    }

    public Expr constant(long value) { return node("const " + value, false, -1); }

    public Expr addVar(long... domain) {
        StringBuilder b = new StringBuilder("V");
        longs(b, domain);
        line(b.toString());
        int var = vars++;
        return node("var " + var, false, var);
    }

    public ListVar addListVar(int owners, int[][] eligible, long... durations) {
        if (durations.length != 0 && durations.length != eligible.length)
            throw new IllegalArgumentException("每个元素需要一个时长");
        StringBuilder b = new StringBuilder("L ").append(owners).append(' ').append(eligible.length);
        for (int i = 0; i < eligible.length; i++) {
            b.append("\nI ").append(durations.length == 0 ? 0 : durations[i]);
            for (int owner : eligible[i]) b.append(' ').append(owner);
        }
        line(b.toString());
        return new ListVar(this, groups++, eligible.length);
    }

    public static final class ListVar {
        private final Model model;
        private final int group, items;

        private ListVar(Model model, int group, int items) {
            this.model = model;
            this.group = group;
            this.items = items;
        }

        public void precede(int before, int after) {
            model.line("P " + group + " " + before + " " + after);
        }

        public Expr load(int owner, long... demands) {
            StringBuilder b = new StringBuilder("load ").append(group).append(' ').append(owner);
            longs(b, demands);
            return model.node(b.toString(), false, -1);
        }

        /** 距离表包含全部元素以及排在最后的仓库。 */
        public Expr distance(int owner, long[][] matrix) {
            if (matrix.length != items + 1) throw new IllegalArgumentException("距离矩阵行数错误");
            StringBuilder b = new StringBuilder("distance ").append(group).append(' ').append(owner);
            for (long[] row : matrix) {
                if (row.length != matrix.length) throw new IllegalArgumentException("距离矩阵列数错误");
                longs(b, row);
            }
            return model.node(b.toString(), false, -1);
        }

        public Expr finish(int item) { return model.node("finish " + group + " " + item, false, -1); }

        /** 指定列表中的 after 是否紧跟 before，两端均为本组实际元素。 */
        public Expr adjacent(int owner, int before, int after) {
            if (before < 0 || before >= items || after < 0 || after >= items)
                throw new IllegalArgumentException("相邻关系的两端必须是本组元素");
            return model.node("adjacent " + group + " " + owner + " " + before + " " + after, true, -1);
        }
    }

    public Expr choose(Expr condition, Expr yes, Expr no) {
        if (!condition.predicate) throw new IllegalArgumentException("choose 需要布尔条件");
        return node("if " + id(condition) + " " + id(yes) + " " + id(no), false, -1);
    }
    public Expr choose(Expr condition, long yes, long no) {
        return choose(condition, constant(yes), constant(no));
    }

    private Expr aggregate(String op, Expr[] terms) {
        StringBuilder b = new StringBuilder(op);
        for (Expr term : terms) b.append(' ').append(id(term));
        return node(b.toString(), false, -1);
    }
    public Expr sum(Expr... terms) { return aggregate("sum", terms); }
    public Expr max(Expr... terms) { return aggregate("max", terms); }
    public Expr count(Expr... conditions) {
        Expr[] terms = new Expr[conditions.length];
        for (int i = 0; i < terms.length; i++) terms[i] = choose(conditions[i], 1, 0);
        return sum(terms);
    }
    public Expr groupCount(Expr[] keys, long key) {
        return groupCount(keys, constant(key));
    }
    public Expr groupCount(Expr[] keys, Expr target) {
        id(target);
        Expr[] conditions = new Expr[keys.length];
        for (int i = 0; i < keys.length; i++) conditions[i] = keys[i].eq(target);
        return count(conditions);
    }

    public Expr groupSum(Expr[] keys, Expr[] values, long key) {
        return groupSum(keys, values, constant(key));
    }
    public Expr groupSum(Expr[] keys, Expr[] values, Expr target) {
        if (keys.length != values.length) throw new IllegalArgumentException("分组键与数值数量不同");
        id(target);
        Expr zero = constant(0);
        Expr[] terms = new Expr[keys.length];
        for (int i = 0; i < keys.length; i++) terms[i] = choose(keys[i].eq(target), values[i], zero);
        return sum(terms);
    }

    public Expr lookup(Expr index, long... table) {
        StringBuilder b = new StringBuilder("lookup ").append(id(index));
        longs(b, table);
        return node(b.toString(), false, -1);
    }

    public void addConstr(Expr condition, String name) { addConstr(condition, name, 1, true); }
    public void addConstr(Expr condition, String name, long weight, boolean hard) {
        addPenalty(choose(condition, 0, 1), name, weight, hard);
    }
    public void addPenalty(Expr violation, String name) { addPenalty(violation, name, 1, false); }
    public void addPenalty(Expr violation, String name, long weight, boolean hard) {
        if (name == null || name.isEmpty() || names.contains(name))
            throw new IllegalArgumentException("规则名必须非空且唯一");
        line("R " + hex(name) + " " + (hard ? "H" : "S") + " " + weight + " " + id(violation));
        names.add(name);
    }
    public void setObjective(Expr e) { check(); id(e); objective = e; }

    public Result optimize() { return optimize(42, 50000, -1, 1, "tabu", "ff"); }
    public synchronized Result optimize(long seed, long steps, long timeMs, int threads,
                                       String algo, String construct) {
        check();
        if (!Arrays.asList("tabu", "sa", "la").contains(algo) ||
                !Arrays.asList("ff", "bf").contains(construct))
            throw new IllegalArgumentException("未知搜索算法或构造方式");
        StringBuilder request = new StringBuilder(lines);
        if (objective != null) request.append("O ").append(objective.node).append('\n');
        request.append("C ").append(Long.toUnsignedString(seed)).append(' ').append(steps).append(' ')
            .append(timeMs).append(' ').append(threads).append(' ').append(algo).append(' ').append(construct);
        byte[] reply = Native.run(request.toString().getBytes(StandardCharsets.US_ASCII));
        String data = new String(reply, 1, reply.length - 1, StandardCharsets.UTF_8);
        if (reply[0] != 0) { failure = data; throw new IllegalArgumentException(data); }
        return new Result(this, data);
    }

    @Override public synchronized void close() {
        closed = true;
        lines.setLength(0);
        names.clear();
    }

    public static final class Result {
        public final long hard, soft, steps;
        public final boolean complete;
        public final String status, reason;
        public final Map<String, Long> penalties;
        private final Model model;
        private final List<Long> values = new ArrayList<Long>();
        private final Map<Long, int[]> sequences = new HashMap<Long, int[]>();

        private static long key(int group, int owner) {
            return ((long)group << 32) | (owner & 0xffffffffL);
        }

        private Result(Model model, String text) {
            this.model = model;
            Map<String, Long> p = new LinkedHashMap<String, Long>();
            String[] summary = null;
            String termination = null;
            for (String line : text.split("\n")) {
                String[] f = line.split(" ");
                switch (f[0]) {
                case "S": summary = f; break;
                case "T": termination = f[1]; break;
                case "V":
                    for (int i = 1; i < f.length; i++) values.add(f[i].equals("_") ? null : Long.valueOf(f[i]));
                    break;
                case "L":
                    int[] seq = new int[f.length - 3];
                    for (int i = 0; i < seq.length; i++) seq[i] = Integer.parseInt(f[i + 3]);
                    sequences.put(key(Integer.parseInt(f[1]), Integer.parseInt(f[2])), seq);
                    break;
                case "R": p.put(unhex(f[1]), Long.valueOf(f[2])); break;
                default: throw new IllegalStateException("未知 Rust 结果行");
                }
            }
            if (summary == null || termination == null) throw new IllegalStateException("求解结果不完整");
            hard = Long.parseLong(summary[1]);
            soft = Long.parseLong(summary[2]);
            steps = Long.parseLong(summary[3]);
            status = summary[4];
            complete = summary[5].equals("1");
            reason = termination;
            penalties = Collections.unmodifiableMap(p);
        }

        public boolean feasible() { return complete && hard == 0; }
        public Long value(Expr var) {
            if (var.model != model || var.var < 0) throw new IllegalArgumentException("需要同一模型的分配变量");
            return values.get(var.var);
        }
        public int[] sequence(ListVar group, int owner) {
            if (group.model != model) throw new IllegalArgumentException("需要同一模型的列表变量");
            int[] seq = sequences.get(key(group.group, owner));
            if (seq == null) throw new IllegalArgumentException("列表不存在");
            return seq.clone();
        }
    }
}
