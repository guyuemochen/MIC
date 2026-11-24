package benchmark;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.IntStream;
import org.optaplanner.core.api.domain.entity.PlanningEntity;
import org.optaplanner.core.api.domain.solution.*;
import org.optaplanner.core.api.domain.valuerange.ValueRangeProvider;
import org.optaplanner.core.api.domain.variable.PlanningListVariable;
import org.optaplanner.core.api.domain.variable.PlanningVariable;
import org.optaplanner.core.api.score.buildin.hardsoftlong.HardSoftLongScore;

/** 输入为双方共享的纯整数事实。输入顺序就是稳定实体编号。 */
public final class Domain {
    public record Data(String kind, int[] header, int[][] rows) {
        static Data read(Path path) throws IOException {
            List<String[]> lines = Files.readAllLines(path).stream().map(String::trim)
                    .filter(s -> !s.isEmpty() && !s.startsWith("#"))
                    .map(s -> s.split("\\s+")).toList();
            int[][] rows = lines.stream().skip(1).map(Data::ints).toArray(int[][]::new);
            return new Data(lines.get(0)[0], ints(lines.get(0)), rows);
        }

        static int[] ints(String[] words) {
            return Arrays.stream(words, 1, words.length).mapToInt(Integer::parseInt).toArray();
        }

        int maxColumn(int col) {
            return Arrays.stream(rows).mapToInt(row -> row[col]).max().orElse(-1);
        }
    }

    public static List<Integer> range(int n) {
        return IntStream.range(0, n).boxed().toList();
    }

    @PlanningSolution
    public abstract static class Plan {
        @ProblemFactProperty public Data data;
        @PlanningScore public HardSoftLongScore score;
    }

    @PlanningEntity
    public static class Shift {
        public int id;
        @ValueRangeProvider(id="qualified") public List<Integer> qualified;
        @PlanningVariable(valueRangeProviderRefs="qualified") public Integer worker;
    }

    @PlanningSolution
    public static class Roster extends Plan {
        @PlanningEntityCollectionProperty public List<Shift> shifts = new ArrayList<>();
    }

    @PlanningEntity
    public static class Lesson {
        public int id;
        @PlanningVariable(valueRangeProviderRefs="times") public Integer time;
        @PlanningVariable(valueRangeProviderRefs="rooms") public Integer room;
    }

    @PlanningSolution
    public static class School extends Plan {
        @ValueRangeProvider(id="times") public List<Integer> times;
        @ValueRangeProvider(id="rooms") public List<Integer> rooms;
        @PlanningEntityCollectionProperty public List<Lesson> lessons = new ArrayList<>();
    }

    @PlanningEntity
    public static class Owner {
        public int id;
        @PlanningListVariable(valueRangeProviderRefs="items") public List<Integer> order = new ArrayList<>();
    }

    @PlanningSolution
    public static class Sequences extends Plan {
        @ValueRangeProvider(id="items") @ProblemFactCollectionProperty public List<Integer> items;
        @PlanningEntityCollectionProperty public List<Owner> owners = new ArrayList<>();
    }

    public static Plan empty(Data d) {
        int n = d.rows.length;
        if (d.kind.equals("ROSTER")) {
            Roster sol = new Roster();
            sol.data = d;
            for (int i = 0; i < n; i++) {
                Shift s = new Shift();
                s.id = i;
                int[] costs = d.rows[i];
                s.qualified = IntStream.range(0, d.header[0]).filter(e -> costs[1+e] >= 0).boxed().toList();
                sol.shifts.add(s);
            }
            return sol;
        }
        if (d.kind.equals("TIMETABLE")) {
            School sol = new School();
            sol.data = d;
            sol.times = range(d.header[0]);
            sol.rooms = range(d.header[1]);
            for (int i = 0; i < n; i++) {
                Lesson lesson = new Lesson();
                lesson.id = i;
                sol.lessons.add(lesson);
            }
            return sol;
        }
        if (!Set.of("DELIVERY","FACTORY").contains(d.kind)) throw new IllegalArgumentException(d.kind);
        Sequences sol = new Sequences();
        sol.data = d;
        sol.items = range(n);
        for (int i = 0; i < d.header[0]; i++) {
            Owner owner = new Owner();
            owner.id = i;
            sol.owners.add(owner);
        }
        return sol;
    }
}
