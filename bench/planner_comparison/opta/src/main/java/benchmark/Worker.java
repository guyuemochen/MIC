package benchmark;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.nio.file.Path;
import java.time.Duration;
import java.util.*;
import benchmark.Domain.*;
import org.optaplanner.core.api.score.director.ScoreDirector;
import org.optaplanner.core.api.solver.*;
import org.optaplanner.core.config.constructionheuristic.ConstructionHeuristicPhaseConfig;
import org.optaplanner.core.config.constructionheuristic.placer.QueuedValuePlacerConfig;
import org.optaplanner.core.config.heuristic.selector.move.composite.UnionMoveSelectorConfig;
import org.optaplanner.core.config.heuristic.selector.move.generic.list.*;
import org.optaplanner.core.config.heuristic.selector.value.ValueSelectorConfig;
import org.optaplanner.core.config.localsearch.LocalSearchPhaseConfig;
import org.optaplanner.core.config.localsearch.decider.acceptor.LocalSearchAcceptorConfig;
import org.optaplanner.core.config.localsearch.decider.forager.LocalSearchForagerConfig;
import org.optaplanner.core.config.score.director.ScoreDirectorFactoryConfig;
import org.optaplanner.core.config.solver.*;
import org.optaplanner.core.impl.heuristic.selector.common.decorator.SelectionFilter;
import org.optaplanner.core.impl.heuristic.selector.move.generic.list.*;

/** 一个进程可以依次接收多条试验请求，以便在相同 JVM 内预热和测量。 */
public final class Worker {
    public static final class MachineFilter implements SelectionFilter<Plan,Object> {
        public boolean accept(ScoreDirector<Plan> sd, Object move) {
            if (move instanceof ListAssignMove<?> a)
                return sd.getWorkingSolution().data.rows()[(Integer)a.getMovedValue()][0]==((Owner)a.getDestinationEntity()).id;
            if (move instanceof ListChangeMove<?> c) return c.getSourceEntity()==c.getDestinationEntity();
            if (move instanceof ListSwapMove<?> s) return s.getLeftEntity()==s.getRightEntity();
            throw new IllegalStateException("未知列表移动 " + move.getClass());
        }
    }

    static Solver<Plan> solver(Plan p, long seed, long ms) {
        Class<?> entity = p instanceof Roster ? Shift.class : p instanceof School ? Lesson.class : Owner.class;
        ConstructionHeuristicPhaseConfig ch = new ConstructionHeuristicPhaseConfig();
        LocalSearchPhaseConfig ls = new LocalSearchPhaseConfig();
        LocalSearchAcceptorConfig accept = new LocalSearchAcceptorConfig();
        if (p instanceof Sequences) accept.withValueTabuSize(8);
        else accept.withEntityTabuSize(8);
        ls.withAcceptorConfig(accept).withForagerConfig(new LocalSearchForagerConfig().withAcceptedCountLimit(32));
        if (p.data.kind().equals("FACTORY")) {
            // 使用官方队列构造器，先按输入顺序取工序，枚举合格机器上的插入位置。
            ch.withEntityPlacerConfig(new QueuedValuePlacerConfig()
                    .withValueSelectorConfig(new ValueSelectorConfig("order").withId("order"))
                    .withMoveSelectorConfig(new ListChangeMoveSelectorConfig()
                            .withValueSelectorConfig(new ValueSelectorConfig().withMimicSelectorRef("order"))
                            .withFilterClass(MachineFilter.class)));
            ls.withMoveSelectorConfig(new UnionMoveSelectorConfig().withMoveSelectors(
                    new ListChangeMoveSelectorConfig().withFilterClass(MachineFilter.class),
                    new ListSwapMoveSelectorConfig().withFilterClass(MachineFilter.class)));
        }
        SolverConfig cfg = new SolverConfig()
                .withSolutionClass(p.getClass()).withEntityClasses(entity)
                .withScoreDirectorFactory(new ScoreDirectorFactoryConfig().withIncrementalScoreCalculatorClass(Score.class))
                .withEnvironmentMode(EnvironmentMode.REPRODUCIBLE).withMoveThreadCount("NONE")
                .withRandomSeed(seed).withTerminationSpentLimit(Duration.ofMillis(ms)).withPhases(ch, ls);
        return SolverFactory.<Plan>create(cfg).buildSolver();
    }

    record Snapshot(long hard,long soft,boolean complete,List<Integer> values,List<List<Integer>> lists) {
        static Snapshot copy(Plan p) {
            List<Integer> values = new ArrayList<>();
            List<List<Integer>> lists = new ArrayList<>();
            if (p instanceof Roster r) r.shifts.forEach(s -> values.add(s.worker));
            if (p instanceof School s) for (Lesson l : s.lessons) {
                values.add(l.time);
                values.add(l.room);
            }
            if (p instanceof Sequences s) s.owners.forEach(o -> lists.add(new ArrayList<>(o.order)));
            boolean complete = p.score != null && p.score.isSolutionInitialized();
            return new Snapshot(p.score==null?0:-p.score.hardScore(),p.score==null?0:-p.score.softScore(),
                    complete, values, lists);
        }
    }

    static final class Observe {
        long start, budget, first = -1;
        Snapshot best, witness;
        void event(Plan plan,long now) {
            long t=now-start;
            if(t>budget) return;
            best=Snapshot.copy(plan);
            if(first<0 && best.complete && best.hard==0) {
                first=t;
                witness=best;
            }
        }
    }

    public static void run(String[] args) throws Exception {
        Plan p = Domain.empty(Data.read(Path.of(args[0])));
        long seed=Long.parseLong(args[1]), ms=Long.parseLong(args[2]);
        Solver<Plan> solver=solver(p,seed,ms);
        Observe obs=new Observe();
        obs.budget=ms*1_000_000;
        obs.best=Snapshot.copy(p);
        solver.addEventListener(e -> obs.event(e.getNewBestSolution(), System.nanoTime()));
        obs.start=System.nanoTime();
        Plan result=solver.solve(p);
        long stop=System.nanoTime();
        if (stop-obs.start<=obs.budget) obs.event(result,stop);
        Snapshot b=obs.best;
        System.out.println("RESULT {\"hard\":"+b.hard+",\"soft\":"+b.soft+",\"complete\":"+b.complete+
                ",\"first_ms\":"+(obs.first<0?"null":Double.toString(obs.first/1e6))+
                ",\"elapsed_ms\":"+((stop-obs.start)/1e6)+",\"values\":"+b.values+",\"lists\":"+b.lists+
                ",\"first_values\":"+(obs.witness==null?"null":obs.witness.values)+
                ",\"first_lists\":"+(obs.witness==null?"null":obs.witness.lists)+"}");
    }

    public static void main(String[] args) throws Exception {
        if (args.length>0) { run(args); return; }
        BufferedReader lines = new BufferedReader(new InputStreamReader(System.in));
        for (String line; (line=lines.readLine())!=null;) {
            if (!line.isBlank()) run(line.trim().split("\\s+"));
        }
    }
}
