import java.util.Arrays;
import org.micplanner.Model;

/** Java 8 真机接口验收，无需测试框架。 */
public final class BindingTest {
    private static void check(boolean value) {
        if (!value) throw new AssertionError("规划器绑定断言失败");
    }

    public static void main(String[] args) {
        Model m = new Model();
        Model.Expr x;
        Model.Result snapshot;
        try (Model autoClose = m) {
            x = m.addVar(0, 1);
            Model.Expr y = m.addVar(0, 1);
            m.addConstr(x.ne(y), "不同时间");
            m.addPenalty(m.count(x.eq(0), y.eq(0)), "夜班负担", 2, false);
            m.setObjective(m.lookup(x, 5, 0));
            snapshot = m.optimize(42, 20, -1, 1, "tabu", "bf");
            check(snapshot.feasible() && snapshot.soft == 2 && snapshot.value(x) == 1);
            check(snapshot.penalties.get("夜班负担") == 2 && snapshot.reason.equals("Steps"));
            m.setObjective(x.mul(5));
            check(m.optimize(42, 20, -1, 1, "tabu", "bf").value(x) == 0);
        }
        check(snapshot.value(x) == 1);
        m.close();
        try { m.optimize(); throw new AssertionError("已关闭模型仍可求解"); }
        catch (IllegalStateException expected) {}
        try (Model timing = new Model()) {
            Model.ListVar g = timing.addListVar(2, new int[][] {{0}, {1}}, 2, 3);
            g.precede(1, 0);
            timing.addConstr(g.load(1, 1, 3).eq(3), "载重");
            timing.setObjective(timing.max(g.finish(0), g.finish(1)).add(
                    g.distance(0, new long[][] {{0,2,1}, {2,0,8}, {3,9,0}})));
            Model.Result result = timing.optimize(42, 10, -1, 1, "tabu", "ff");
            check(result.feasible() && result.soft == 9 && Arrays.equals(result.sequence(g, 0), new int[]{0}));
        }
        try (Model invalid = new Model()) {
            invalid.addVar();
            try { invalid.optimize(); throw new AssertionError("空值域未报错"); }
            catch (IllegalArgumentException expected) {}
            try { invalid.optimize(); throw new AssertionError("错误状态未保留"); }
            catch (IllegalStateException expected) {}
        }
        try (Model grouping = new Model()) {
            Model.Expr[] keys = {grouping.addVar(0), grouping.addVar(1), grouping.addVar(0)};
            Model.Expr[] values = {grouping.constant(7), grouping.constant(11), grouping.constant(13)};
            grouping.addConstr(grouping.groupCount(keys, 0).eq(2), "人数");
            grouping.setObjective(grouping.groupSum(keys, values, 0));
            check(grouping.optimize(42, 0, -1, 1, "tabu", "ff").soft == 20);
        }
        try (Model routing = new Model()) {
            Model.Expr x2 = routing.addVar(2);
            Model.ListVar route = routing.addListVar(2, new int[][] {{0}, {0}});
            route.precede(0, 1);
            Model.Expr forward = route.adjacent(0, 0, 1);
            Model.Expr reverse = route.adjacent(0, 1, 0);
            Model.Expr empty = route.adjacent(1, 0, 1);
            routing.addConstr(forward, "顺序");
            routing.addConstr(reverse, "逆序", 3, false);
            routing.setObjective(routing.choose(empty, x2.mul(4), x2)
                    .add(routing.choose(forward, x2, routing.constant(0))));
            try { route.adjacent(0, 0, 2); throw new AssertionError("仓库被当作相邻元素"); }
            catch (IllegalArgumentException expected) {}
            Model.Result r = routing.optimize(42, 10, -1, 1, "tabu", "ff");
            check(Arrays.equals(r.sequence(route, 0), new int[]{0,1}) && r.sequence(route, 1).length == 0);
            check(r.hard == 0 && r.soft == 7 && r.penalties.get("顺序") == 0 && r.penalties.get("逆序") == 3);
        }
    }
}
