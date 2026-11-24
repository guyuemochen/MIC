package benchmark;

import java.util.Arrays;
import java.util.List;
import org.optaplanner.core.api.score.buildin.hardsoftlong.HardSoftLongScore;
import org.optaplanner.core.api.score.calculator.IncrementalScoreCalculator;
import benchmark.Domain.*;

/** OptaPlanner 原生增量计分。正数业务代价在框架中取负，按硬分、软分排序。 */
public final class Score implements IncrementalScoreCalculator<Plan, HardSoftLongScore> {
    private Data d;
    private long hard, soft;
    private int[] load;
    private int[][] workerAt, teacherAt, cohortAt, roomAt;
    private JobTimes jobs;

    public void resetWorkingSolution(Plan plan) {
        d = plan.data;
        hard = soft = 0;
        jobs = null;
        if (plan instanceof Roster r) {
            load = new int[d.header()[0]];
            workerAt = new int[d.header()[0]][d.maxColumn(0)+1];
            r.shifts.forEach(s -> scalar(s, 1));
        } else if (plan instanceof School s) {
            int periods = d.header()[0];
            teacherAt = new int[d.maxColumn(0)+1][periods];
            cohortAt = new int[d.maxColumn(1)+1][periods];
            roomAt = new int[d.header()[1]][periods];
            s.lessons.forEach(l -> scalar(l, 1));
        } else {
            Sequences seq = (Sequences) plan;
            if (d.kind().equals("FACTORY")) jobs = new JobTimes(d);
            load = new int[seq.owners.size()];
            for (Owner o : seq.owners) afterListVariableChanged(o, "order", 0, o.order.size());
        }
    }

    private void pair(int[] counts, int k, int sign) {
        // k 个已占用实体再加一个，恰多 k 对冲突。
        counts[k] += sign;
        hard += sign > 0 ? counts[k] - 1 : -counts[k];
    }

    private void capacity(int owner, int delta) {
        int cap = d.header()[1];
        hard -= Math.max(0, load[owner] - cap);
        load[owner] += delta;
        hard += Math.max(0, load[owner] - cap);
    }

    private void scalar(Object entity, int sign) {
        if (entity instanceof Shift shift && shift.worker != null) {
            int worker = shift.worker;
            int[] row = d.rows()[shift.id];
            capacity(worker, sign);
            pair(workerAt[worker], row[0], sign);
            soft += (long) sign * row[worker + 1];
        } else if (entity instanceof Lesson lesson && lesson.time != null) {
            int t = lesson.time;
            int[] row = d.rows()[lesson.id];
            pair(teacherAt[row[0]], t, sign);
            pair(cohortAt[row[1]], t, sign);
            if (lesson.room != null) pair(roomAt[lesson.room], t, sign);
            soft += (long) sign * Math.abs(t-row[2]);
        }
    }

    public void beforeVariableChanged(Object entity, String variable) { scalar(entity, -1); }
    public void afterVariableChanged(Object entity, String variable) { scalar(entity, 1); }
    public void beforeEntityAdded(Object entity) {}
    public void afterEntityAdded(Object entity) {}
    public void beforeEntityRemoved(Object entity) {}
    public void afterEntityRemoved(Object entity) {}

    private int distance(int left, int right) {
        int x1 = left < 0 ? d.header()[2] : d.rows()[left][0];
        int y1 = left < 0 ? d.header()[3] : d.rows()[left][1];
        int x2 = right < 0 ? d.header()[2] : d.rows()[right][0];
        int y2 = right < 0 ? d.header()[3] : d.rows()[right][1];
        return Math.abs(x1-x2) + Math.abs(y1-y2);
    }

    private void route(Owner owner, int from, int to, int sign) {
        List<Integer> items = owner.order;
        int prev = from == 0 ? -1 : items.get(from-1);
        int demand = 0;
        long edges = 0;
        for (int i = from; i < to; i++) {
            int item = items.get(i);
            edges += distance(prev, item);
            demand += d.rows()[item][2];
            prev = item;
        }
        edges += distance(prev, to == items.size() ? -1 : items.get(to));
        soft += sign * edges;
        capacity(owner.id, sign * demand);
    }

    public void beforeListVariableChanged(Object entity, String name, int from, int to) {
        Owner owner = (Owner) entity;
        if (jobs == null) route(owner, from, to, -1);
        else jobs.detach(owner.order, from, to);
    }

    public void afterListVariableChanged(Object entity, String name, int from, int to) {
        Owner owner = (Owner) entity;
        if (jobs == null) route(owner, from, to, 1);
        else jobs.attach(owner.order, from, to);
    }

    public HardSoftLongScore calculateScore() {
        if (jobs != null) return jobs.score();
        return HardSoftLongScore.of(-hard, -soft);
    }

    /** 变动区段和原有紧后工序为起点，仅把新的最长路径传播到相应后继。 */
    private static final class JobTimes {
        private final Data data;
        private final int n;
        private final int[][] jobNext;
        private final int[] prev, next, degree, queue;
        private final boolean[] present, dirty, affected;
        private final long[] end;
        private boolean hadCycle;

        JobTimes(Data data) {
            this.data = data;
            n = data.rows().length;
            prev = new int[n]; next = new int[n]; degree = new int[n]; queue = new int[n];
            Arrays.fill(prev, -1);
            Arrays.fill(next, -1);
            end = new long[n];
            present = new boolean[n]; dirty = new boolean[n]; affected = new boolean[n];
            int[] counts = new int[n];
            for (int[] row : data.rows()) if (row[2] >= 0) counts[row[2]]++;
            jobNext = new int[n][];
            for (int i=0; i<n; i++) jobNext[i] = new int[counts[i]];
            Arrays.fill(counts, 0);
            for (int i=0; i<n; i++) {
                int p = data.rows()[i][2];
                if (p>=0) jobNext[p][counts[p]++] = i;
            }
        }

        void detach(List<Integer> list, int from, int to) {
            if (from>0) next[list.get(from-1)] = -1;
            for (int i=from; i<to; i++) {
                int v = list.get(i);
                prev[v] = next[v] = -1;
                present[v] = false; dirty[v] = true;
            }
            if (to<list.size()) {
                int v = list.get(to);
                prev[v] = -1; dirty[v] = true;
            }
        }

        void attach(List<Integer> list, int from, int to) {
            int left = from == 0 ? -1 : list.get(from-1);
            for (int i=from; i<to; i++) {
                int v = list.get(i);
                prev[v] = left;
                if (left>=0) next[left]=v;
                present[v]=true; dirty[v]=true; left=v;
            }
            int right = to==list.size() ? -1 : list.get(to);
            if (left>=0) next[left]=right;
            if (right>=0) {prev[right]=left; dirty[right]=true;}
        }

        HardSoftLongScore score() {
            int count=0;
            for(int i=0;i<n;i++) {
                affected[i]=dirty[i] || hadCycle;
                dirty[i]=false;
                if (affected[i]) queue[count++]=i;
            }
            for (int k=0;k<count;k++) {
                int v=queue[k];
                for (int s: jobNext[v]) if (present[s] && !affected[s]) {affected[s]=true; queue[count++]=s;}
                int s=next[v];
                if (s>=0 && !affected[s]) {affected[s]=true; queue[count++]=s;}
            }
            int active=0, tail=0;
            for(int i=0;i<n;i++) {
                if (!affected[i]) continue;
                if (!present[i]) {end[i]=0; continue;}
                active++;
                int p=data.rows()[i][2], m=prev[i];
                degree[i]=(p>=0 && present[p] && affected[p]?1:0)+(m>=0 && affected[m]?1:0);
                if (degree[i]==0) queue[tail++]=i;
            }
            for(int head=0;head<tail;head++) {
                int v=queue[head], p=data.rows()[v][2], m=prev[v];
                end[v]=data.rows()[v][1]+Math.max(p<0 || !present[p]?0:end[p],m<0?0:end[m]);
                for(int s:jobNext[v]) if(present[s] && affected[s] && --degree[s]==0) queue[tail++]=s;
                int s=next[v];
                if(s>=0 && affected[s] && --degree[s]==0) queue[tail++]=s;
            }
            hadCycle=tail!=active;
            if(hadCycle) return HardSoftLongScore.of(-1,0);
            long longest=0;
            for(long t:end) longest=Math.max(longest,t);
            return HardSoftLongScore.of(0,-longest);
        }
    }
}
