import java.util.StringJoiner;
import org.micplanner.Model;
import org.micplanner.Model.Expr;
import org.micplanner.Model.ListVar;
import org.micplanner.Model.Result;

/** 分配、排课、路线与加工顺序的 Java 8 建模示例。 */
public final class PlannerCases {
    private static void solveAndPrint(Model model, String name, Expr[] variables,
                                      ListVar group, int owners) {
        Result result = model.optimize(7, 120, -1, 1, "tabu", "bf");
        StringJoiner values = new StringJoiner(",");
        for (Expr v : variables) values.add(String.valueOf(result.value(v)));
        StringJoiner routes = new StringJoiner("/");
        for (int owner = 0; owner < owners; owner++) {
            StringJoiner route = new StringJoiner(",");
            for (int item : result.sequence(group, owner)) route.add(String.valueOf(item));
            routes.add(route.toString());
        }
        StringJoiner penalties = new StringJoiner(",");
        for (Long penalty : result.penalties.values()) penalties.add(penalty.toString());
        System.out.println(String.join("\t", name, Long.toString(result.hard),
                Long.toString(result.soft), Long.toString(result.steps), result.status,
                result.reason, result.complete ? "1" : "0", values.toString(),
                routes.toString(), penalties.toString()));
    }

    private static void roster() {
        try (Model model = new Model()) {
            Expr[] shifts = {model.addVar(0, 1), model.addVar(0, 1), model.addVar(1)};
            for (int employee = 0; employee < 2; employee++) {
                Expr[] assigned = new Expr[shifts.length];
                for (int i = 0; i < shifts.length; i++) assigned[i] = shifts[i].eq(employee);
                model.addConstr(model.count(assigned).le(2), "capacity" + employee);
            }
            for (int i = 0; i < shifts.length; i++)
                model.addConstr(shifts[i].eq(1), "pref" + i, 2, false);
            solveAndPrint(model, "roster", shifts, null, 0);
        }
    }

    private static void timetable() {
        try (Model model = new Model()) {
            Expr[] time = new Expr[3], room = new Expr[3], variables = new Expr[6];
            for (int i = 0; i < time.length; i++) {
                time[i] = variables[2 * i] = model.addVar(0, 1);
                room[i] = variables[2 * i + 1] = model.addVar(0, 1);
            }
            model.addConstr(time[0].ne(time[1]), "teacher");
            for (int i = 0; i < time.length; i++) {
                for (int j = i + 1; j < time.length; j++) {
                    Expr conflict = time[i].eq(time[j]).mul(room[i].eq(room[j]));
                    model.addPenalty(conflict, "room" + i + j, 1, true);
                }
            }
            model.setObjective(model.sum(time));
            solveAndPrint(model, "timetable", variables, null, 0);
        }
    }

    private static void delivery() {
        try (Model model = new Model()) {
            ListVar trucks = model.addListVar(2, new int[][] {{0, 1}, {0, 1}, {0, 1}});
            int[] locations = {1, 3, 8, 0};  // 最后一个位置是仓库。
            long[][] matrix = new long[locations.length][locations.length];
            for (int i = 0; i < locations.length; i++)
                for (int j = 0; j < locations.length; j++)
                    matrix[i][j] = Math.abs(locations[i] - locations[j]);
            for (int v = 0; v < 2; v++)
                model.addConstr(trucks.load(v, 1, 1, 1).le(2), "capacity" + v);
            model.setObjective(model.sum(trucks.distance(0, matrix), trucks.distance(1, matrix)));
            solveAndPrint(model, "delivery", new Expr[0], trucks, 2);
        }
    }

    private static void factory() {
        try (Model model = new Model()) {
            ListVar machines = model.addListVar(2, new int[][] {{0}, {1}, {1}, {0}}, 2, 1, 2, 1);
            machines.precede(0, 1);
            machines.precede(2, 3);
            Expr[] finish = new Expr[4];
            for (int i = 0; i < finish.length; i++) finish[i] = machines.finish(i);
            model.setObjective(model.max(finish));
            solveAndPrint(model, "factory", new Expr[0], machines, 2);
        }
    }

    public static void main(String[] args) {
        roster();
        timetable();
        delivery();
        factory();
    }
}
